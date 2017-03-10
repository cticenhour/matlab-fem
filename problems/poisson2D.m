% Finite element solution of Poisson's equation (-U_xx - U_yy = 1) on a
% unit square with dirichlet boundary conditions (U_at_boundary = 0)
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% Original: September 2016
% Revised: March 2017

%=============================

MOOSE_comparison = 0;   % requires output data CSV file from MOOSE
                        % requires MOOSE and MATLAB codes have same mesh!

%=============================
    
% Generate 2D mesh using DistMesh (see 
% http://persson.berkeley.edu/distmesh for more info)
% KEY NOTES FOR DISTMESH
% node_list = node positions -- N-by-2 array contains x,y coordinates
%                               for each of the N nodes
% triangle_list = triangle indices -- row associated with each triangle 
%                                     has 3 integer entries to specify 
%                                     node numbers corresponding to 
%                                     rows in node_list

len = 1;
wid = 1;
initial_edge_width = min(wid,len)/10;

geo_dist_func = @(p) drectangle(p,0,wid,0,len);

bounds = [0,0;wid,len];
important_pts = [0,0;wid,0;0,len;wid,len];

% Note: @huniform refers to uniform mesh
% See http://persson.berkeley.edu/distmesh/ for usage info
[node_list, triangle_list] = distmesh2d(geo_dist_func,@huniform,...
    initial_edge_width,bounds,important_pts);

boundary_edges = boundedges(node_list,triangle_list);
edge_nodes = unique(boundary_edges);

num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);
    
% Initialize parts of system KU=F, where u is solution vector
K = zeros(num_nodes,num_nodes);
F = zeros(num_nodes,1);

% Weak form components
K = buildKlaplacian(K,triangle_list,node_list,edge_nodes,1);
F = buildFsource(F,triangle_list,node_list,edge_nodes,1);
% Boundary conditions
[K,F] = pecBC(K,F,edge_nodes);

% Solve system of equations
U = K\F;

% Plot solution
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U,...
    'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar
if vectorized_method == 1
    title('FEM solution, vectorized')
elseif vectorized_method == 0
    title('FEM solution with loops')
end
xlabel('X')
ylabel('Y')

% Generate and plot analytic solution on triangular mesh
U_analytic = analyticFXN_poisson2D(node_list,50);

figure
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_analytic,'edgecolor','k','facecolor','interp')
view(2),axis equal,colorbar
title('Analytic solution to 50 terms in infinite series');
xlabel('X')
ylabel('Y')

RMS_analytic_error = sqrt(sum((U_analytic - U).^2)/length(U))
%two_norm_error = norm(U_analytic-U,2)

% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_analytic-U,'edgecolor','k','facecolor','interp')
% view(2),axis equal,colorbar
% title('Absolute error')

if MOOSE_comparison == 1
    sortMOOSEoutput
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_MOOSE,'edgecolor','k','facecolor','interp')
    view(2),axis equal,colorbar
    title('MOOSE solution');
    xlabel('X')
    ylabel('Y')
    
    RMS_MOOSE_error = sqrt(sum((U_MOOSE - U).^2)/length(U))
end