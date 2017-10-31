% Finite element solution of the scalar wave equation with two components
% in a waveguide with perfectly conducting walls
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% October 2017
% last update: October 31, 2017


%=============================
% PROBLEM PARAMETERS AND CONTANTS
%=============================

% load physical constants
phys_constants

waveguide_width = 10;
waveguide_length = 80;

num_components = 3;
dim = 2;

m = 1;

init_E = 1;

f = 20e6;
omega = 2*pi*f;

k0 = omega/c;

k = k0*(1+1i*0);

beta = sqrt(k^2 - (pi/waveguide_width)^2);

source = 0;

%=============================
%   MESH IMPORT
%=============================

mesh = 'C:/Users/cticenhour/Code/matlab-fem/meshes/waveguide.msh';

% Create usable connectivity info from GMSH mesh file

[node_list,triangle_list,boundary_edges,boundary_names] = ...
                                    gmsh2matlab2d(mesh);
% Create edge nodes and edges arrays, labeled with appropriate names
total_bounds = length(boundary_names);
for i = 1:total_bounds
    eval([boundary_names{i},'_edge_nodes = nonzeros(unique(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,1:2)));'])
    eval([boundary_names{i},'_edges = [nonzeros(ismember(boundary_edges(:,3),',...
        num2str(i),').*boundary_edges(:,1)), nonzeros(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,2))];'])
end

edge_nodes = unique(boundary_edges);

% Clean up task - Remove duplicate corners on walls from exit and port
port_edge_nodes(ismember(port_edge_nodes,intersect(top_edge_nodes,port_edge_nodes))) = [];
port_edge_nodes(ismember(port_edge_nodes,intersect(top_edge_nodes,exit_edge_nodes))) = [];

exit_edge_nodes(ismember(exit_edge_nodes,intersect(bottom_edge_nodes,port_edge_nodes))) = [];
exit_edge_nodes(ismember(exit_edge_nodes,intersect(bottom_edge_nodes,exit_edge_nodes))) = [];

% Create other needed things
wall_edges = [top_edges; bottom_edges];
wall_edge_nodes = [top_edge_nodes; bottom_edge_nodes];
num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);

%=============================
%   BUILD SYSTEM
%=============================

% Initialize parts of system KU=F, where U is solution vecor for all three
% components
K = zeros(num_nodes*3,num_nodes*3);
F = zeros(num_nodes*3,1);

Kx = zeros(num_nodes,num_nodes);
Fx = zeros(num_nodes,1);

Ky = zeros(num_nodes,num_nodes);
Fy = zeros(num_nodes,1);

Kz = zeros(num_nodes,num_nodes);
Fz= zeros(num_nodes,1);

% Weak Form X
% GradDiv
K = buildKgraddiv(K,triangle_list,node_list,0,0,0);
% Laplacian
Kx = buildKlaplacian(Kx,triangle_list,node_list,0);
% Coefficient * Field
Kx = buildKcoeff(Kx,triangle_list,node_list,0,k0);
% Right hand source term
Fx = buildFsource(Fx,triangle_list,node_list,0,source);

% Weak Form Y
% Laplacian
Ky = buildKlaplacian(Ky,triangle_list,node_list,0);
% Coefficient * Field
Ky = buildKcoeff(Ky,triangle_list,node_list,0,k0);
% Right hand source term
Fy = buildFsource(Fy,triangle_list,node_list,0,source);

% Weak Form Z
% Laplacian
Kz = buildKlaplacian(Kz,triangle_list,node_list,0);
% Coefficient * Field
Kz = buildKcoeff(Kz,triangle_list,node_list,0,k0);
% Right hand source term
Fz = buildFsource(Fz,triangle_list,node_list,0,source);

% Boundary Conditions X

% Boundary Conditions Y 
% Port BC at z = 0
[Ky,Fy] = portBC(Ky,Fy,triangle_list,node_list,port_edges, port_edge_nodes,k0,waveguide_width,1);
% Absorbing BC at z = 80
[Ky,Fy] = absorbingBC(Ky,Fy,triangle_list,node_list,exit_edges,exit_edge_nodes,k0,waveguide_width);
% Natural (E' = 0) condition on walls
[Ky,Fy] = neumannBC(Ky,Fy,triangle_list,node_list,wall_edges,wall_edge_nodes,k0,waveguide_width,0);

% Boundary Conditions Z 
% Port BC at z = 0
[Kz,Fz] = portBC(Kz,Fz,triangle_list,node_list,port_edges, port_edge_nodes,k0,waveguide_width,2);
% Absorbing BC at z = 80
[Kz,Fz] = absorbingBC(Kz,Fz,triangle_list,node_list,exit_edges,exit_edge_nodes,k0,waveguide_width);
% PEC (E = 0) condition on walls
[Kz,Fz] = pecBC(Kz,Fz,wall_edge_nodes);

K(1:num_nodes,1:num_nodes) = Ky; 
K((num_nodes+1):(num_nodes*2),(num_nodes+1):(num_nodes*2)) = Kz;

F(1:num_nodes,1) = Fy;
F((num_nodes+1):(num_nodes*2),1) = Fz;

% Solve system of equations
U = K\F;

% extract components from solution vector
Ey = U(1:num_nodes,1);
Ez = U((num_nodes+1):(num_nodes*2),1);


figure
subplot(3,1,1)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),abs(Ez),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Magnitude of E_z')
xlabel('Z (m)')
ylabel('Y (m)')

% solution real component
subplot(3,1,2)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(Ez),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Real part of E_z')
xlabel('Z (m)')
ylabel('Y (m)')

% solution imaginary component
subplot(3,1,3)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),imag(Ez),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Imaginary part of E_z')
xlabel('Z (m)')
ylabel('Y (m)')
        

figure
subplot(3,1,1)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),abs(Ey),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Magnitude of E_y')
xlabel('Z (m)')
ylabel('Y (m)')

% solution real component
subplot(3,1,2)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(Ey),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Real part of E_y')
xlabel('Z (m)')
ylabel('Y (m)')

% solution imaginary component
subplot(3,1,3)
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),imag(Ey),...
    'edgecolor','none','facecolor','interp');
view(2)
axis image
colorbar
title('Imaginary part of E_y')
xlabel('Z (m)')
ylabel('Y (m)')

figure 
spy(K)
