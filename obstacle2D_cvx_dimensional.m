% OBSTACLE2D_CVX_DIMENSIONAL
% 
% Two-dimensional cross-sectional model of Institute Ice Stream dynamics
% using convex optimization (CVX) with dimensional physical units.
%
% This script solves the ice flow problem using a variational formulation
% that minimizes the total power dissipation subject to constraints.
% The model includes:
%   - Glen's flow law (n=3, p=4/3) for ice rheology
%   - Basal sliding with spatially variable friction
%   - Temperature-dependent viscosity (thermal-mechanical coupling)
%   - Finite element discretization on unstructured triangular mesh
%
% Dependencies:
%   - CVX optimization toolbox (http://cvxr.com/cvx/)
%   - distmesh2d for mesh generation (http://persson.berkeley.edu/distmesh/)
%   - vel_profile_full.mat (observation data)
%   - iceColorMap.mat (colormap for visualization, optional)
%
% Output:
%   - Velocity field u (m/s)
%   - Temperature field T (degrees C)
%   - Visualization plots comparing model to observations
%
% Reference: See accompanying paper for full methodology and results

clc
clear

%% ========================================================================
% INITIALIZATION
% ========================================================================
% Physical parameters
p = 4/3;                    % Flow law exponent (Glen's law: n=3, p=4/3)
rho = 900;                  % Ice density [kg/m^3]
g = 9.8;                    % Gravitational acceleration [m/s^2]
alpha = 2.4E-3;             % Surface slope [rad]
f = rho*g*sin(alpha);       % Driving stress [Pa/m]

% Domain parameters
dx = 880;                   % Nominal grid spacing [m]
vert_scale = 0.375;         % Vertical scaling factor for domain
dy = vert_scale*dx;         % Vertical grid spacing [m]
% Basal friction parameters (tau_c)
% These define the spatially variable basal friction law
% Multiple numerical experiments are commented out below - uncomment to use
%
% Numerical experiment #1 (uniform strength)
% beta = 4.1E6;
% rock_trans = 15.1E3;
% sed_const = 23.95E3;
% sed_var = 0E3;
% sed_trans_l = -23E3;
% sed_trans_r = 72E3;
% ch_str = 0;
% ch_decay = 7E3;
% ch_loc = 0E3;
%
% Numerical experiment #2 (linear strength)
% beta = 4E6;
% rock_trans = 13.5E3;
% sed_const = 16.55E3;
% sed_var = 16E3;
% sed_trans_l = -23E3;
% sed_trans_r = 71E3;
% ch_str = 0;
% ch_decay = 7E3;
% ch_loc = 0E3;
%
% Numerical experiment #3 (all plastic)
% beta = 0;
% rock_trans = -23E3;
% sed_const = 17E3;
% sed_var = 15E3;
% sed_trans_l = -23E3;
% sed_trans_r = 71E3;
% ch_str = 500;
% ch_decay = 20E3;
% ch_loc = 13E3;

% Numerical experiment #4 (rock/sediment with channel) - ACTIVE
beta_var = 5E6;             % Variation in rock friction coefficient [Pa^(5/2) s^(1/2) m^(-5/2)]
beta = 2.9E6;                % Base rock friction coefficient [Pa^(5/2) s^(1/2) m^(-5/2)]
rock_trans = 10E3;           % Transition from rock to sediment [m]
sed_const = 19.4E3;          % Constant sediment friction [Pa s/m]
sed_var = 8E3;               % Variable sediment friction [Pa s/m]
sed_trans_l = -23E3;         % Left sediment transition [m]
sed_trans_r = 71E3;          % Right sediment transition [m]
ch_str = 24.5E3;             % Channel strength [Pa s/m]
ch_decay = 5.9E3;            % Channel decay length scale [m]
ch_loc = 10E3;               % Channel location [m]

% Basal friction law: tau_c(x,y,u)
% Piecewise function defining friction as function of position and velocity:
%   - Left boundary (x < sed_trans_l): very high friction (1E8 Pa s/m)
%   - Rock region (sed_trans_l < x < rock_trans): power law friction
%   - Sediment region (rock_trans < x < sed_trans_r): linear + channel friction
%   - Right boundary (x > sed_trans_r): very high friction
%   - Only active below bed surface (y < 450 + dy/2)
tau_c =@(x,y,u) (...
                 heaviside(sed_trans_l - x).*1E8.*abs(u) + ...
                 heaviside(x - sed_trans_l).*heaviside(rock_trans - x).*(beta_var*(x/80E3) + beta).^2.*pow_abs(u,5/2) + ...
                 heaviside(x - rock_trans).*(sed_var*((80E3-x)/80E3) + ...
                                             ch_str*(exp(-abs(x-ch_loc)/ch_decay)) + ...
                                             sed_const + ...
                                             1E8*heaviside(x - sed_trans_r)).*abs(u)... 
                 ).*heaviside((450+dy/2) - y);


%% ========================================================================
% MESH GENERATION
% ========================================================================
% Define domain polygon vertices (scaled by 40 km)
% Domain represents cross-section of Institute Ice Stream
pv = 40E3.*[0,1;2,1;2,.93;.3,.9;0,.925;-.75,.925;-.75,1;0,1];

% Generate unstructured triangular mesh using distmesh2d
% @dpoly: polygon distance function
% @huniform: uniform mesh size function
% dx: target mesh spacing
% Bounding box: 40E3.*[-.85,.75;2.1,1.1]
[xy,t] = distmesh2d(@dpoly,@huniform,dx,40E3.*[-.85,.75;2.1,1.1],pv,pv);

% Apply vertical scaling to mesh coordinates
xy(:,2) = vert_scale*((xy(:,2)-40E3) - min(xy(:,2)-40E3));

% Mesh statistics
nN = size(xy,1);             % Number of nodes
nE = size(t,1);              % Number of elements (triangles)
b = unique(boundedges(xy,t)); % Boundary node indices
nB = size(b,1);              % Number of boundary nodes

% Identify boundary edges (edges that appear only once)
e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];  % All edges
e = sort(e,2);                           % Sort edge endpoints
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:);     % Boundary edges (appear once)


%% ========================================================================
% TEMPERATURE INITIALIZATION AND ITERATIVE THERMAL-MECHANICAL COUPLING
% ========================================================================
% Initialize temperature field with linear profile from surface to bed
% Surface temperature: 248 K (-25 C)
% Temperature gradient: 25 K / 1500 m
T = 25/1500*(1500-xy(:,2)) + 248;
T_old = T;
res = 1;                     % Residual for convergence check

% Iterative thermal-mechanical coupling loop
while(res > 1E-3)
    % Under-relaxation for temperature update
    T = .3*T + (1-.3)*T_old;
    
    %% Temperature-dependent ice rheology
    % Glen's flow law: epsilon_dot = A * tau^n
    % where A depends on temperature via Arrhenius relation
    
    A = 3.5e-25;              % Preexponential constant [1/s*Pa^3]
    E = 1;                    % Enhancement factor (dimensionless)
    Q_h = 0*xy(:,1);          % Melting point adjusted activation energy [J/mol]
    R = 8.314;                % Gas constant [J/K*mol]
    Tstar = 263.15;           % Activation threshold temperature [K] (-10 C)
    p0 = 7e-8;                % Pressure heating coefficient [K/Pa]
    
    % Pressure-dependent temperature (melting point depression)
    P = rho*g*(1500-xy(:,2)); % Hydrostatic pressure [Pa]
    T_h = T + p0*P;           % Pressure-adjusted temperature [K]
    Tstar_h = Tstar*(0*xy(:,1) + 1) + p0*P;  % Pressure-adjusted threshold [K]
    
    % Activation energy depends on temperature regime
    Q2 = 60e3;                % Lower activation energy [J/mol] (T < Tstar)
    Q3 = 115e3;               % Higher activation energy [J/mol] (T >= Tstar)
    Q_h(heaviside(T_h-Tstar_h)==0) = Q2;
    Q_h(heaviside(T_h-Tstar_h)==1) = Q3;
    
    % Compute temperature-dependent viscosity coefficient
    % a = (1/2) * A^(-1/n) where n=3, so a = (1/2) * A^(-1/3)
    a = 1/2*(A*E*exp((-Q_h./R).*((1./T_h)-(1./Tstar_h)))).^(-1/3);

    %% ====================================================================
    % BUILD FINITE ELEMENT SYSTEM
    % ====================================================================
    % Initialize sparse matrices for finite element discretization
    % A, B: gradient operators (x and y components)
    % D: element-to-node interpolation matrix
    % F: force vector (driving stress)
    A = sparse(nE,nN);
    B = sparse(nE,nN);
    D = zeros(nE,nN);
    F = zeros(1,nN);
    
    tau_area = zeros(nE,1);   % Element areas for integration
    
    % Loop over elements to build system matrices
    for E = 1:nE
        nodes = t(E,:);         % Node indices for this element
        xyE = [ones(3,1),xy(nodes,:)];  % Element coordinate matrix
        Area = abs(det(xyE))/2;         % Element area
        tau_area(E) = Area;
        
        % Compute gradient operators via inverse coordinate transformation
        C = inv(xyE);
        A_E = C(2,:);           % x-gradient operator
        B_E = C(3,:);           % y-gradient operator
        
        % Element shape function integrals (linear elements)
        D_E = 1/3*ones(1,3);    % Constant interpolation
        F_E = f*Area/3*ones(1,3);  % Force vector contribution
        
        % Assemble into global matrices
        A(E,nodes) = A(E,nodes) + A_E;
        B(E,nodes) = B(E,nodes) + B_E;
        D(E,nodes) = D(E,nodes) + D_E;
        F(nodes) = F(nodes) + F_E;
    end
    
    % Compute boundary integration weights (b_dx)
    % These account for boundary element lengths for surface integral
    b_dx = zeros(size(b,1),1);
    for N = 1:size(b,1)
        edges = eB(sum(eB == b(N),2) == 1,:);  % Edges connected to boundary node
        for i = 1:size(edges,1)
            edge = edges(i,:);
            % Add half-length of each connected edge (trapezoidal rule)
            b_dx(N) = b_dx(N) + .5*sqrt(sum(diff(xy(edge,:),1).^2,2));
        end
    end


    %% ====================================================================
    % SOLVE VELOCITY FIELD
    % ====================================================================
    % Minimize total power dissipation using CVX
    % Objective function:
    %   J = (1/p) * int( a * |grad u|^p ) - int( f * u ) + int_bdry( tau_c * |u| )
    %   where p = 4/3 (Glen's law with n=3)
    
    cvx_begin
    cvx_quiet true
        variables u(nN)         % Velocity field at nodes
        
        % Power dissipation functional
        % First term: internal viscous dissipation
        % Second term: work done by driving stress
        % Third term: basal friction dissipation
        obj = (1/p)*sum((D*a).*tau_area.*pow_pos(norms([A*u,B*u],2,2),p)) ...
              - F*u ...
              + sum(b_dx.*tau_c(xy(b,1),xy(b,2),u(b)));
        
        subject to
            u >= 0;             % Non-negative velocity (downslope flow)
            u(xy(:,1) > 80E3 - dx/2) == 0;   % Right boundary: no slip
            u(xy(:,1) < -30E3 + dx/2) == 0;   % Left boundary: no slip
        minimize(obj)
    cvx_end
    
    % Compute velocity gradients and derived quantities
    u_x = A*u;                  % x-component of velocity gradient
    u_y = B*u;                  % y-component of velocity gradient
    
    % Effective viscosity (from Glen's law)
    mu = (D*a).*(sqrt((u_x).^2 + (u_y).^2)).^(1/2);
    
    % Effective stress and strain rate
    tau_E = sqrt((mu.*u_x).^2 + (mu.*u_y).^2);  % Effective stress [Pa]
    epsilon_E = sqrt((u_x/2).^2 + (u_y/2).^2); % Effective strain rate [1/s]
    
    % Frictional heating rate [W/m^3]
    f_therm = 2*tau_E.*epsilon_E;
    
    %% Temperature-dependent thermal conductivity
    k1 = 9.828;                 % Conductivity preexponential [W/m*K]
    k2 = 5.7;                   % Conductivity postexponential [1/K]
    k = k1*exp(-k2*1e-3.*T);    % Thermal conductivity [W/m*K]
    
    % Basal temperature boundary condition
    tau_T =@(y) (55./(1500-y.*heaviside((450+dy/2) - y))).*heaviside((450+dy/2) - y);
    
    % Assemble thermal forcing from frictional heating
    F_therm = zeros(1,nN);
    for E = 1:nE
        nodes = t(E,:);
        xyE = [ones(3,1),xy(nodes,:)];
        Area = abs(det(xyE))/2;
        F_therm_E = f_therm(E)*Area/3*ones(1,3);
        F_therm(nodes) = F_therm(nodes) + F_therm_E;
    end
    
    %% ====================================================================
    % SOLVE TEMPERATURE FIELD
    % ====================================================================
    % Solve heat equation with frictional heating source
    % Minimize: int( k * |grad T|^2 ) - int( f_therm * T ) - int_bdry( k * tau_T * T )
    
    T_old = T;
    cvx_begin
    cvx_quiet true
        variables T(nN)         % Temperature field at nodes
        
        % Heat conduction functional
        obj = sum((D*k).*tau_area.*pow_pos(norms([A*T,B*T],2,2),2)) ...
              - 2E10*F_therm*T ...  % Frictional heating source
              - sum(b_dx.*k(b).*tau_T(xy(b,2)).*T(b));  % Basal boundary condition
        
        subject to
            T <= 273;            % Temperature cannot exceed melting point [K]
            T(xy(:,2) > 1500-dx/10) == 248;  % Surface temperature boundary condition [K]
        minimize(obj)
    cvx_end
    
    % Check convergence
    res = norm(T-T_old)/norm(T);
    disp(res)
end
%% ========================================================================
% VISUALIZATION
% ========================================================================
% Convert units for display
u = u*pi*1E7;                % Convert to m/s (scaling factor)
T = T - 273;                 % Convert to degrees Celsius

% Figure 1: Velocity field
trisurf(t,xy(:,1),xy(:,2),u,u,...
       'edgecolor','none','facecolor','interp');
hold on
% Mark boundary conditions
plot3(xy(b(u(b) < 1e-6),1),xy(b(u(b) < 1e-6),2),max(u)+0*b(u(b) < 1e-6),'r.','MarkerSize',20)
plot3(xy(b(b_dx==0),1),xy(b(b_dx==0),2),max(u)+0*b(b_dx==0),'b.','MarkerSize',20)
view(2), colorbar
title('Velocity Field [m/s]')
xlabel('Distance [m]')
ylabel('Elevation [m]')

% Figure 2: Comparison with observations
load vel_profile_full.mat
figure('Position',[10 10 500 400])

subplot(2,1,1)
plot(profile_path-30.5E3,profile_cross,'LineWidth',3)
hold on
plot(xy(xy(:,2) > 1500-dx/10,1),u(xy(:,2) > 1500-dx/10),'LineWidth',3)
axis([-30E3,80E3,0,400])
xlabel('Distance [m]')
ylabel('Velocity [m/s]')
title('Surface Velocity Profile')
legend('Observed','Modeled','Location','best')

subplot(2,1,2)
plot(profile_path-30.5E3,gradient(profile_cross),'LineWidth',3)
hold on
plot(xy(xy(:,2) > 1500-dx/10,1),gradient(u(xy(:,2) > 1500-dx/10)),'LineWidth',3)
axis([-30E3,80E3,-80,50])
xlabel('Distance [m]')
ylabel('Velocity Gradient [m/s/m]')
title('Velocity Gradient')
legend('Observed','Modeled','Location','best')
set(gca,'FontSize',16)

% Figure 3: Temperature field
figure('Position',[10 10 1000 200])
% Load colormap if available, otherwise use default
if exist('iceColorMap.mat','file')
    load iceColorMap
    colormap(iceColorMap)
else
    colormap('parula')  % Use MATLAB default colormap if iceColorMap not available
end
trisurf(t,xy(:,1),xy(:,2),T,T,...
       'edgecolor','none','facecolor','interp');
axis off
view(2)
colorbar
title('Temperature Field [°C]')

% Figure 4: Basal friction parameterization
figure('Position',[10 10 500 150])
yyaxis left
plot(sed_trans_l:1E2:rock_trans,beta+beta_var*(sed_trans_l:1E2:rock_trans)/80E3,'LineWidth',3)
axis([-30E3,80E3,1E6,4E6])
ylabel('Rock Friction Coefficient')
yyaxis right
plot(rock_trans:1E2:sed_trans_r,sed_var*((80E3-(rock_trans:1E2:sed_trans_r))/80E3) + ch_str*(exp(-abs((rock_trans:1E2:sed_trans_r)-ch_loc)/ch_decay)) + sed_const,'LineWidth',3)
axis([-30E3,80E3,15E3,55E3])
xlabel('Distance [m]')
ylabel('Sediment Friction [Pa s/m]')
title('Basal Friction Parameterization')
set(gca,'FontSize',16)

%figure 
%bedmap2_profile(profile_lat(profile_path < 110E3),profile_lon(profile_path < 110E3))