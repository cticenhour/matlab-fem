% Finite element solution of a simple wave equation (curl(curl(E)) + E = 1)
% in a waveguide with perfectly conducting walls (E_wall = 0)
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% December 2016
% last update: December 6, 2016

clear all
close all

%=============================
% SWITCHES
%=============================

mesh_program = 0; % distMesh = 1, Gmsh = 0
    % NOTE: gmsh2matlab REQUIRES NODE, ELEMENT, AND EDGE TXT FILES TO BE 
    %       CREATED FROM GMSH OUTPUT FILE

MOOSE_comparison = 0;   % requires output data CSV file from MOOSE

%=============================

if mesh_program == 1
    
    % Generate 2D mesh using DistMesh (see 
    % http://persson.berkeley.edu/distmesh for more info)
    % KEY NOTES FOR DISTMESH
    % node_list = node positions -- N-by-2 array contains x,y coordinates
    %                               for each of the N nodes
    % triangle_list = triangle indices -- row associated with each triangle 
    %                                     has 3 integer entries to specify 
    %                                     node numbers corresponding to 
    %                                     rows in node_list

    len = 5;
    wid = 1;
    initial_edge_width = min(len,wid)/10;

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
elseif mesh_program == 0
    gmsh2matlab_waveguide
end
    
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
    
    centroid_x = (1/3)*(node_list(current_nodes(1),1) + ...
        node_list(current_nodes(2),1) + node_list(current_nodes(3),1));
    
    centroid_y = (1/3)*(node_list(current_nodes(1),2) + ...
        node_list(current_nodes(2),2) + node_list(current_nodes(3),2));
        
    % CONSTRUCT K WITHOUT BOUNDARY CONDITIONS

    for r = 1:3
        for s = r:3
            % Determine if nodes r and s are interior or on boundary
            node_r = current_nodes(r);
            node_s = current_nodes(s);
            is_r_on_boundary = sum(top_edge_nodes == node_r) + ...
                sum(bottom_edge_nodes == node_r) + ...
                sum(side_edge_nodes == node_r);
            is_s_on_boundary = sum(top_edge_nodes == node_s) + ...
                sum(bottom_edge_nodes == node_s) + ...
                sum(side_edge_nodes == node_s);

            % Determine gradients in XY from coefficients of basis fxn
            gradX_r = linear_basis_coefficients(2,r);
            gradY_r = linear_basis_coefficients(3,r);
            gradX_s = linear_basis_coefficients(2,s);
            gradY_s = linear_basis_coefficients(3,s);

            laplacian = triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);
            
            % Estimate second term using one point gaussian quadrature
            trial_r = linear_basis_coefficients(1,r)+...
                linear_basis_coefficients(2,r)*centroid_x + ...
                linear_basis_coefficients(3,r)*centroid_y;
            
            trial_s = linear_basis_coefficients(1,s)+...
                linear_basis_coefficients(2,s)*centroid_x + ...
                linear_basis_coefficients(3,s)*centroid_y;
            
            linear = (triangle_area/3)*trial_r*trial_s;
            
            integral = laplacian + linear;

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
        is_rF_on_boundary = sum(side_edge_nodes == node_rF) + ...
            sum(top_edge_nodes == node_rF) + ...
            sum(bottom_edge_nodes == node_rF);
        
        if is_rF_on_boundary == 0

            increment = triangle_area/3;
            F(node_rF,1) = F(node_rF,1) + increment;
        end
    end
end  

% Add boundary conditions to K and F

% Set sides BC
for j = 1:length(side_edge_nodes)
    BC_current_node = side_edge_nodes(j);
    K(BC_current_node,BC_current_node) = 1;
    F(BC_current_node,1) = 0;
end

% Set bottom BC

for j = 1:length(bottom_edge_nodes)
    BC_current_node = bottom_edge_nodes(j);
    K(BC_current_node,BC_current_node) = 1;
    F(BC_current_node,1) = 1;
end

% Set top BC

for j = 1:length(top_edge_nodes)
    BC_current_node = top_edge_nodes(j);
    K(BC_current_node,BC_current_node) = 1;
    F(BC_current_node,1) = 1;
end

% Solve system of equations
U = K\F;

% Plot solution
trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U,...
    'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar
xlabel('X')
ylabel('Y')

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
        