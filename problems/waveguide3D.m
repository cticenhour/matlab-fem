% Finite element solution of the scalar wave equation with three components
% in a waveguide with perfectly conducting walls
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% September 2017
% last update: September 20, 2017
    
mesh = '/home/cticenhour/projects/EELS/meshes/waveguide.msh';

% load physical constants
run(constants)

%=============================
% IMPORTANT PROBLEM PARAMETERS
%=============================

waveguide_width = 10;
waveguide_length = 80;

m = 1;

init_E = 1;

omega = 2*pi*20e6;

k0 = omega/c;

k = k0*(1+1i*0);

beta = sqrt(k^2 - (pi/waveguide_width)^2);

source = 0;

%=============================
%   MESH IMPORT
%=============================

% Create usable connectivity info from GMSH mesh file

[node_list,triangle_list,boundary_edges,boundary_names] = ...
                                    gmsh2matlab2d(mesh);
% Create edge nodes and edges arrays, labeled with appropriate names
total_bounds = length(boundary_names);
for i = 1:total_bounds
    eval([boundary_names{i},'_edge_nodes = nonzeros(unique(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,1:2)));'])
    eval([boundary_names{i},'_edges = [nonzeros(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,1)), nonzeros(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,2))];'])
end

edge_nodes = unique(boundary_edges);

% Clean up task - Remove duplicate corners on walls from exit and port
port_edge_nodes(ismember(port_edge_nodes,intersect(top_edge_nodes,port_edge_nodes))) = [];
port_edge_nodes(ismember(port_edge_nodes,intersect(top_edge_nodes,exit_edge_nodes))) = [];

exit_edge_nodes(ismember(exit_edge_nodes,intersect(bottom_edge_nodes,port_edge_nodes))) = [];
exit_edge_nodes(ismember(exit_edge_nodes,intersect(bottom_edge_nodes,exit_edge_nodes))) = [];

% Create other needed things
wall_edge_nodes = [top_edge_nodes; bottom_edge_nodes];
num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);

%=============================
%   BUILD SYSTEM
%=============================

% Initialize parts of system KU=F, where U is solution vector
K = zeros(num_nodes,num_nodes);
F = zeros(num_nodes,1);

% Weak Form 
% Term 1 - Laplacian (E_xx + E_yy)
K = buildKlaplacian(K,triangle_list,node_list,wall_edge_nodes);
% Term 2 - Coefficient * Field
K = buildKcoeff(K,triangle_list,node_list,wall_edge_nodes,k0);
% Right hand source term
F = buildFsource(F,triangle_list,node_list,wall_edge_nodes,source);

% Boundary Conditions
% Port BC at z = 0
[K,F] = portBC(K,F,triangle_list,node_list,port_edges, port_edge_nodes,k0,waveguide_width);
% Absorbing BC at z = 80
[K,F] = absorbingBC(K,F,triangle_list,node_list,exit_edges,exit_edge_nodes,k0,waveguide_width);
% PEC (E = 0) condition on walls
[K,F] = pecBC(K,F,wall_edge_nodes);

% Solve system of equations
U = K\F;


        