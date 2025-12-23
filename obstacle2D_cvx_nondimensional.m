% OBSTACLE2D_CVX_NONDIMENSIONAL
% 
% Two-dimensional cross-sectional model of Institute Ice Stream dynamics
% using convex optimization (CVX) with nondimensional units.
%
% This script solves the ice flow problem using a variational formulation
% that minimizes the total power dissipation subject to constraints.
% The model includes:
%   - Glen's flow law (n=3, p=4/3) for ice rheology
%   - Basal sliding with spatially variable friction
%   - Finite element discretization on unstructured triangular mesh
%
% This nondimensional version uses scaled units and is useful for:
%   - Understanding scaling relationships
%   - Testing numerical methods
%   - Comparing with analytical solutions
%
% Dependencies:
%   - CVX optimization toolbox (http://cvxr.com/cvx/)
%   - distmesh2d for mesh generation (http://persson.berkeley.edu/distmesh/)
%   - vel_profiles.mat (observation data for comparison)
%
% Output:
%   - Velocity field u (nondimensional)
%   - Visualization plots comparing model to observations
%
% Reference: See accompanying paper for full methodology and results

clc
clear

%% ========================================================================
% INITIALIZATION
% ========================================================================
% Grid spacing (nondimensional)
dx = 0.025;                  % Nominal grid spacing

% Physical parameters (nondimensional)
f = 1.;                      % Driving stress (normalized to unity)

% Basal friction law: tau_c(x,y,u)
% Piecewise function defining friction as function of position and velocity:
%   - Left region (x < 0.4): power law friction with coefficient 133
%   - Right region (x > 0.4): linear friction with step function
%   - Only active below bed surface (y < 1 - dx/2)
tau_c =@(x,y,u) (...
                 heaviside(.4 - x).*133.*pow_abs(u,5/2) + ...
                 heaviside(x - .4).*(2E-5 + 1.1*heaviside(x - 1.8)).*abs(u)...
                ).*heaviside((1-dx/2) - y);


%% ========================================================================
% MESH GENERATION
% ========================================================================
% Define domain polygon vertices (nondimensional coordinates)
% Domain represents simplified cross-section
pv = [0,1;2,1;2,.93;.3,.9;0,.93;0,1];

% Generate unstructured triangular mesh using distmesh2d
% @dpoly: polygon distance function
% @huniform: uniform mesh size function
% dx: target mesh spacing
% Bounding box: [-.1,.75;2.1,1.1]
[xy,t] = distmesh2d(@dpoly,@huniform,dx,[-.1,.75;2.1,1.1],pv,pv);

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
% BUILD FINITE ELEMENT SYSTEM
% ========================================================================
% Initialize sparse matrices for finite element discretization
% A, B: gradient operators (x and y components)
% D: element-to-node interpolation matrix
% F: force vector (driving stress)
A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(nE,nN);
F = zeros(1,nN);

tau_area = zeros(nE,1);     % Element areas for integration

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
    D_E = Area/3*ones(1,3); % Area-weighted interpolation
    F_E = Area/3*ones(1,3); % Force vector contribution
    
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


%% ========================================================================
% SOLVE VELOCITY FIELD
% ========================================================================
% Minimize total power dissipation using CVX
% Objective function:
%   J = (3/4) * int( |grad u|^(4/3) ) - int( f * u ) + int_bdry( tau_c * |u| )
%   where p = 4/3 (Glen's law with n=3)
%   Note: viscosity coefficient absorbed into scaling

cvx_begin
    variables u(nN)         % Velocity field at nodes
    
    % Power dissipation functional
    % First term: internal viscous dissipation (nondimensional)
    % Second term: work done by driving stress
    % Third term: basal friction dissipation
    obj = 3/4*sum(tau_area.*pow_pos(norms([A*u,B*u],2,2),4/3)) ...
          - F*u ...
          + sum(b_dx.*tau_c(xy(b,1),xy(b,2),u(b)));
    
    subject to
        u >= 0;             % Non-negative velocity (downslope flow)
    minimize(obj)
cvx_end


%% ========================================================================
% VISUALIZATION
% ========================================================================
% Figure 1: Velocity field
trisurf(t,xy(:,1),xy(:,2),u,u,...
       'edgecolor','none','facecolor','interp');
hold on
% Mark boundary conditions
plot3(xy(b(u(b) < 1e-6),1),xy(b(u(b) < 1e-6),2),max(u)+0*b(u(b) < 1e-6),'r.','MarkerSize',20)
plot3(xy(b(b_dx==0),1),xy(b(b_dx==0),2),max(u)+0*b(b_dx==0),'b.','MarkerSize',20)
view(2), axis equal, colorbar
title('Velocity Field (Nondimensional)')
xlabel('Distance (nondimensional)')
ylabel('Elevation (nondimensional)')

% Figure 2: Comparison with observations
figure
% Scale velocity for comparison (scaling factor 350)
plot(xy(xy(:,2) > 1-dx/2,1)./2,350*u(xy(:,2) > 1-dx/2)/max(u(xy(:,2) > 1-dx/2)),'LineWidth',3)
hold on
load vel_profiles.mat
j = 6;  % Profile index
% Scale and shift observation profile for comparison
plot(.12+.825*profile_path(:,j)./max(profile_path(:,j)),profile_cross(:,j),'LineWidth',3)
xlabel('Distance (scaled)')
ylabel('Velocity (scaled)')
title('Surface Velocity Profile Comparison')
legend('Modeled','Observed','Location','best')

% Figure 3: Velocity gradient comparison
figure
plot(xy(xy(:,2) > 1-dx/2,1)./2,gradient(u(xy(:,2) > 1-dx/2)/max(u(xy(:,2) > 1-dx/2))))
hold on
plot(.12+.825*profile_path(:,j)./max(profile_path(:,j)),gradient(profile_cross(:,j)/max(profile_cross(:,j))),'LineWidth',3)
xlabel('Distance (scaled)')
ylabel('Velocity Gradient (normalized)')
title('Velocity Gradient Comparison')
legend('Modeled','Observed','Location','best')
