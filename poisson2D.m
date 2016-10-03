% Finite element solution of Poisson's equation (-U_xx - U_yy = 1) on a
% unit square with dirichlet boundary conditions (U_at_boundary = 0)
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% September 2016

clear all
close all

% Switch to determine vectorized or non-vectorized solution method

vectorized_method = 0;

% Generate 2D mesh using DistMesh (see http://persson.berkeley.edu/distmesh/
% for more info)
% KEY NOTES FOR DISTMESH
% node_list = node positions -- N-by-2 array contains x,y coordinates for 
%                               each of the N nodes
% triangle_list = triangle indices -- row associated with each triangle has 
%                                     3 integer entries to specify node 
%                                     numbers corresponding to rows in 
%                                     node_list

len = 1;
wid = 1;
initial_edge_width = wid/10;

geo_dist_func = @(p) drectangle(p,0,wid,0,len);

bounds = [0,0;len,wid];
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

% Applying element-oriented FEM algorithm to construct K and F
for i = 1:num_triangles
    
    current_nodes = triangle_list(i,:);
    
    % FOR LINEAR LAGRANGE BASIS FUNCTIONS (\phi = a + bx + cy)
    triangle_coord_matrix = [1, node_list(current_nodes(1,1),:);...
                             1, node_list(current_nodes(1,2),:);...
                             1, node_list(current_nodes(1,3),:);];
                         
    linear_basis_coefficients = inv(triangle_coord_matrix);
    triangle_area = abs(det(triangle_coord_matrix))/2;
    
    if vectorized_method == 1
        
        gradients = linear_basis_coefficients(2:3,:);
        
        K_element = triangle_area*(gradients')*gradients;
        
        F_element = triangle_area/3;
        
        K(current_nodes,current_nodes) = K(current_nodes,current_nodes)...
            + K_element;
        F(current_nodes) = F(current_nodes) + F_element;
        
    elseif vectorized_method == 0
        
        % CONSTRUCT K WITHOUT BOUNDARY CONDITIONS
    
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

                    if r ~= s
                        K(node_s,node_r) = K(node_s,node_r) + integral;   
                    end

                elseif is_r_on_boundary == 1 && is_s_on_boundary == 0

                    K(node_s,node_r) = K(node_s,node_r) + integral;   

                elseif is_r_on_boundary == 0 && is_s_on_boundary == 1

                    K(node_r,node_s) = K(node_r,node_s) + integral;

                end
            end
        end
        
        % CONSTRUCT F WITHOUT BOUNDARY CONDITIONS
    
        for r = 1:3
            node_rF = current_nodes(r);
            is_rF_on_boundary = sum(edge_nodes == node_rF);

            if is_rF_on_boundary == 0

                increment = triangle_area/3;
                F(node_rF,1) = F(node_rF,1) + increment;
            end
        end
    end  
end

% Add boundary conditions to K and F

if vectorized_method == 1
    
    K(edge_nodes,:) = 0;     % Zero boundary rows of K
    F(edge_nodes) = 0;       % Zero boundary rows of F   
    K(edge_nodes,edge_nodes) = eye(length(edge_nodes),length(edge_nodes));
    % Set K such that U at all boundary nodes has the condition U = 0
    
elseif vectorized_method == 0
    
    for j = 1:length(edge_nodes)
        BC_current_node = edge_nodes(j);
        K(BC_current_node,BC_current_node) = 1;
        F(BC_current_node,1) = 0;
    end
    
end

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

RMS_error = sqrt(sum((U_analytic - U).^2)/length(U))
%two_norm_error = norm(U_analytic-U,2)

% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_analytic-U,'edgecolor','k','facecolor','interp')
% view(2),axis equal,colorbar
% title('Absolute error')

        