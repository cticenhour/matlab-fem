% Finite element solution of -U_xx - U_yy = 1 on a unit square with
% dirichlet boundary conditions U_boundary = 0
% VERSION 2 - generalized away from MATLAB specific routines and functions
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
% Laboratory
% September 27, 2016

clear all
close all

% Use DistMesh to create 2D triangular mesh
length = 1;
width = 1;
initial_edge_width = width/10;

geo_dist_func = @(p) drectangle(p,0,width,0,length);

bounds = [0,0;length,width];
important_pts = [0,0;width,0;0,length;width,length];

% Note: @huniform refers to uniform mesh
% See http://persson.berkeley.edu/distmesh/ for usage info
[node_list, triangle_list] = distmesh2d(geo_dist_func,@huniform,...
    initial_edge_width,bounds,important_pts);

boundary_edges = boundedges(node_list,triangle_list);
edge_nodes = unique(boundary_edges);

num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);

% Initialize parts of system Ku=f, where u is solution vector
K = zeros(num_nodes,num_nodes);
f = zeros(num_nodes,1);

% Applying element-oriented FEM algorithm to construct K and f 
for i = 1:num_triangles
    
    current_nodes = triangle_list(i,:);
    
    % FOR LINEAR LAGRANGE BASIS FUNCTIONS (\phi = a + bx + cy)
    triangle_coord_matrix = [1, node_list(current_nodes(1,1),:);...
                             1, node_list(current_nodes(1,2),:);...
                             1, node_list(current_nodes(1,3),:);];
                         
    linear_basis_coefficients = inv(triangle_coord_matrix);
    triangle_area = abs(det(triangle_coord_matrix))/2;
    
    for r = 1:3
        for s = r:3
            % Determine if nodes r and s are interior or on boundary
            node_r = current_nodes(r);
            node_s = current_nodes(s);
            is_r_on_boundary = sum(edge_nodes == node_r);
            is_s_on_boundary = sum(edge_nodes == node_s);
            
            % Determine gradients in XY from coefficients of basis fxn
            gradX_r = linear_basis_coefficients(2,r);
            gradY_r = linear_basis_coefficients(3,r);
            gradX_s = linear_basis_coefficients(2,s);
            gradY_s = linear_basis_coefficients(3,s);
            
            integral = triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);

            if is_r_on_boundary == 0 && is_s_on_boundary == 0
                
                K(node_r,node_s) = K(node_r,node_s) + integral;
                  
                f(node_r) = f(node_r) + triangle_area/3;
                
                if r ~= s
                    K(node_s,node_r) = K(node_s,node_r) + integral;
                      
                    f(node_s) = f(node_s) + triangle_area/3;
                end
            elseif is_r_on_boundary == 1 && is_s_on_boundary == 0
                K(node_s,node_r) = K(node_s,node_r) + integral;   
                f(node_s) = f(node_s) + triangle_area/3;
                
                K(node_r,node_r) = 1;    % Set K and f such that u(r,1) = 0
                f(node_r) = 0;
            elseif is_r_on_boundary == 0 && is_s_on_boundary == 1
                K(node_r,node_s) = K(node_r,node_s) + integral;
                f(node_r) = f(node_r) + triangle_area/3;
                
                K(node_s,node_s) = 1;   % Set K and f such that u(s,1) = 0
                f(node_s) = 0;
            else
                K(node_r,node_r) = 1;    % Set K and f such that u(r,1) = 0
                f(node_r) = 0;         % and u(s,1) = 0
                K(node_s,node_s) = 1;    
                f(node_s) = 0;
            end
        end
    end
end

% Solve system of equations 
u = K\f;

% Plot solution
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),u,...
    'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar
title('FEM solution')

% Generate and plot analytic solution on triangular mesh
u_analytic = analyticFXN_poisson2D(node_list,50);

figure
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),u_analytic,'edgecolor','k','facecolor','interp')
view(2),axis equal,colorbar
title('Analytic solution to 50 terms in infinite series');

two_norm_error = norm(u_analytic-u,2)

% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),u_analytic-u,'edgecolor','k','facecolor','interp')
% view(2),axis equal,colorbar
% title('Absolute error')