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

mesh_program = 1; % distMesh = 1, Gmsh = 0
    % NOTE: gmsh2matlab REQUIRES NODE, ELEMENT, AND EDGE TXT FILES TO BE 
    %       CREATED FROM GMSH OUTPUT FILE

MOOSE_comparison = 0;   % requires output data CSV file from MOOSE

%=============================
% IMPORTANT CONSTANTS
%=============================

width = 10;
len = 80;

m = 1;

init_E = 1;

source = 0;

omega = 2*pi*20e6;

mu0 = 4*pi*1e-7;

c = 3e8;

Z0 = mu0*c;

k = omega/c*(1+1i*0.1);

current_term = 1i*k*Z0*source;

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

    initial_edge_width = min(len,width)/5;

    geo_dist_func = @(p) drectangle(p,0,width,0,len);

    bounds = [0,0;width,len];
    important_pts = [0,0;width,0;0,len;width,len];

    % Note: @huniform refers to uniform mesh
    % See http://persson.berkeley.edu/distmesh/ for usage info
    [node_list, triangle_list] = distmesh2d(geo_dist_func,@huniform,...
        initial_edge_width,bounds,important_pts);

    boundary_edges = boundedges(node_list,triangle_list);
    edge_nodes = unique(boundary_edges);
   
    edge_coords = node_list(edge_nodes(:),:);
    
    bottom_logic = ismember(edge_coords(:,2),0);
    top_logic = ismember(edge_coords(:,2),len);
    left_side_logic = ismember(edge_coords(:,1),0);
    right_side_logic = ismember(edge_coords(:,1),width);
    
    top_edge_nodes = nonzeros(top_logic.*edge_nodes);
    bottom_edge_nodes = nonzeros(bottom_logic.*edge_nodes);
    left_edge_nodes = nonzeros(left_side_logic.*edge_nodes);
    right_edge_nodes = nonzeros(right_side_logic.*edge_nodes);
    
    % set duplicate side corners to null values, removing them from left and right 
    left_edge_nodes(1) = [];
    left_edge_nodes(length(left_edge_nodes)) = [];
    
    right_edge_nodes(1) = [];
    right_edge_nodes(length(right_edge_nodes)) = [];
    
    side_edge_nodes = [left_edge_nodes; right_edge_nodes];
    
    
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
        
    % CONSTRUCT K CONTRIBUTION FOR CURRENT TRIANGLE WITHOUT BOUNDARY 
    % CONDITIONS

    for r = 1:3
        for s = r:3
            % Determine if nodes r and s are on boundary, neglecting check 
            % for top (side = PEC, bottom = PORT, top = natural)
            node_r = current_nodes(r);
            node_s = current_nodes(s);
            is_r_on_port_or_sides = sum(bottom_edge_nodes == node_r) + ...
                sum(side_edge_nodes == node_r);
            is_r_on_top = sum(top_edge_nodes == node_r);
            is_s_on_port_or_sides = sum(bottom_edge_nodes == node_s) + ...
                sum(side_edge_nodes == node_s);
            is_s_on_top = sum(top_edge_nodes == node_s);

            % Determine gradients in XY from coefficients of basis fxn
            const_r = linear_basis_coefficients(1,r);
            gradX_r = linear_basis_coefficients(2,r);
            gradY_r = linear_basis_coefficients(3,r);
            
            const_s = linear_basis_coefficients(1,s);
            gradX_s = linear_basis_coefficients(2,s);
            gradY_s = linear_basis_coefficients(3,s);

            % Estimate first term using one point gaussian quadrature
            laplacian = triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);
            
            % Estimate second term using one point gaussian quadrature
            trial_r = const_r + gradX_r*centroid_x + gradY_r*centroid_y;
            trial_s = const_s + gradX_s*centroid_x + gradY_s*centroid_y;
            
            linear = triangle_area*(-k^2)*trial_r*trial_s;
            
            LHS = laplacian + linear;

            if is_r_on_port_or_sides == 0 && is_s_on_port_or_sides == 0

                 K(node_r,node_s) = K(node_r,node_s) + LHS;

                if r ~= s
                    K(node_s,node_r) = K(node_s,node_r) + LHS;   
                end

            elseif is_r_on_port_or_sides == 1 && is_s_on_port_or_sides == 0

                K(node_s,node_r) = K(node_s,node_r) + LHS;              

            elseif is_r_on_port_or_sides == 0 && is_s_on_port_or_sides == 1

                K(node_r,node_s) = K(node_r,node_s) + LHS;

            end
            
            % Add absorbing boundary condition to top boundary
            % (first-order)
            
            absorbing_BC_r = triangle_area*(gradY_r*trial_s + 1i*k*trial_r*trial_s);
            
            absorbing_BC_s = triangle_area*(gradY_s*trial_r + 1i*k*trial_s*trial_r);
            
            if is_r_on_top == 1 && is_s_on_top == 0
               
                K(node_r,node_s) = K(node_r, node_s) + absorbing_BC_r;
                
            elseif is_r_on_top == 0 && is_s_on_top == 1
                
                K(node_s,node_r) = K(node_s, node_r) + absorbing_BC_s;
                
            end
        end
    end

    % CONSTRUCT F CONTRIBUTION WITH PORT BOUNDARY CONDITION ON BOTTOM

    for r = 1:3
        node_rF = current_nodes(r);
        is_rF_on_boundary = sum(side_edge_nodes == node_rF) + ...
            sum(top_edge_nodes == node_rF) + ...
            sum(bottom_edge_nodes == node_rF);
        
        is_rF_on_bottom_boundary = sum(bottom_edge_nodes == node_rF);
        
        % Estimate integral via three point gaussian quadrature (incrementing 
        % one node at a time) 
        increment = triangle_area/3*current_term;
        
        if is_rF_on_bottom_boundary == 1 % PORT BOUNDARY CONDITION 
            
            trial_rF = linear_basis_coefficients(1,r)+...
                linear_basis_coefficients(2,r)*centroid_x + ...
                linear_basis_coefficients(3,r)*centroid_y;
            
            F(node_rF,1) = F(node_rF,1) + increment + 1i*k*init_E*sin(pi*m*node_list(node_rF,1)/width)*trial_rF;
            
            if K(node_rF,node_rF) == 0 
                K(node_rF,node_rF) = 1;
            end
        end
        
        if is_rF_on_boundary == 0
            F(node_rF,1) = F(node_rF,1); %+ increment;
        end
    end
    
    % Add side PEC boundary conditions

    for j = 1:length(side_edge_nodes)
        BC_current_node = side_edge_nodes(j);
        K(BC_current_node,BC_current_node) = 1;
        F(BC_current_node,1) = 0;
    end

    % Add side periodic boundary conditions (IN PROGRESS - MIGHT NEED
    % CONTINUITY IN VALUE, BELOW, and CONTINUITY IN THE FIRST DERIVATIVE IN
    % MAIN LOOP)

%     for j = 1:length(right_edge_nodes)
%         BC_current_node_right = right_edge_nodes(j);
%         BC_current_node_left = left_edge_nodes(j);
%         % Continuity in magnitude
%         K(BC_current_node_right,BC_current_node_right) = 1;
%         K(BC_current_node_right,BC_current_node_left) = -1;
%         F(BC_current_node_right,1) = 0;
%     end
    
end  

% Solve system of equations
U = K\F;

phase = atan2(imag(U),real(U));

rF_cycles = 20;
t_max = rF_cycles/omega;

t = 0:t_max/100:t_max;

E_time = U*exp(1i*omega*t);

% Plot solution
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),abs(U),...
%     'edgecolor','k','facecolor','interp');
% view(2),axis equal,colorbar
% title('Magnitude of E_z')
% xlabel('X')
% ylabel('Y')
% 
% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(U),...
%     'edgecolor','k','facecolor','interp');
% view(2),axis equal,colorbar
% title('Real part of E_z')
% xlabel('X')
% ylabel('Y')
% 
% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),imag(U),...
%     'edgecolor','k','facecolor','interp');
% view(2),axis equal,colorbar
% title('Imaginary part of E_z')
% xlabel('X')
% ylabel('Y')
% 
% figure
% trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),phase,...
%     'edgecolor','k','facecolor','interp');
% view(2),axis equal,colorbar
% title('Phase of E_z')
% xlabel('X')
% ylabel('Y')

% Plot time solution
figure
for j=1:length(t)
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(E_time(:,j)),...
    'edgecolor','k','facecolor','interp');
    view(2),axis equal
    h = colorbar;
    set(h,'ylim',[-0.5 0.5])
    xlabel('X')
    ylabel('Y')
    F = getframe(gcf);
end 

% Take a slice from a trisurf plot at the center of the waveguide 
% (x = width/2)
% num_pts = 30;
% method = 'natural';
% xy = [ones(num_pts,1).*(width/2), linspace(0,len,num_pts)'];
% z = griddata(node_list(:,1),node_list(:,2),phase,xy(:,1),xy(:,2));
% 
% figure
% plot(xy(:,2),z);

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
        