% Finite element solution of the scalar wave equation
% in a waveguide with perfectly conducting walls
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% December 2016
% last update: February 2, 2017

clear all
close all

%=============================
% SWITCHES AND PLOTTING OPTIONS
%=============================

mesh_program = 1; % distMesh = 1, Gmsh = 0
    % NOTE: gmsh2matlab REQUIRES NODE, ELEMENT, AND EDGE TXT FILES TO BE 
    %       CREATED FROM GMSH OUTPUT FILE

MOOSE_comparison = 0;   % requires output data CSV file from MOOSE

% Plotting switches
real_E = 0;
imag_E = 0;
mag_E = 1;
phase_E = 0;
time_E = 0;
analytic_time_E = 0;
slice = 0;

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

eps0 = 8.854e-12; % F/m

c = 3e8;

Z0 = mu0*c;

k0 = omega/c;

k = k0*(1+1i*0);

beta = sqrt(k^2 - (pi/width)^2);

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

    geo_dist_func = @(p) drectangle(p,0,len,0,width);

    bounds = [0,0;len,width];
    important_pts = [0,0;len,0;0,width;len,width];

    % Note: @huniform refers to uniform mesh
    % See http://persson.berkeley.edu/distmesh/ for usage info
    [node_list, triangle_list] = distmesh2d(geo_dist_func,@huniform,...
        initial_edge_width,bounds,important_pts);

    % close distmesh rendering of mesh
    close(1)
    
    boundary_edges = boundedges(node_list,triangle_list);
    edge_nodes = unique(boundary_edges);
   
    edge_coords = node_list(edge_nodes(:),:);
    
    port_logic = ismember(edge_coords(:,1),0);
    exit_logic = ismember(edge_coords(:,1),len);
    top_side_logic = ismember(edge_coords(:,2),width);
    bottom_side_logic = ismember(edge_coords(:,2),0);
    
    exit_edge_nodes = nonzeros(exit_logic.*edge_nodes);
    port_edge_nodes = nonzeros(port_logic.*edge_nodes);
    top_edge_nodes = nonzeros(top_side_logic.*edge_nodes);
    bottom_edge_nodes = nonzeros(bottom_side_logic.*edge_nodes);
    
    % set duplicate side corners to null values, removing them from top and
    % bottom
    top_edge_nodes(1) = [];
    top_edge_nodes(length(top_edge_nodes)) = [];
    
    bottom_edge_nodes(1) = [];
    bottom_edge_nodes(length(bottom_edge_nodes)) = [];
    
    wall_edge_nodes = [top_edge_nodes; bottom_edge_nodes];
    
    
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
    
    centroid_x = (1/3)*(node_list(current_nodes(1,1),1) + ...
        node_list(current_nodes(1,2),1) + node_list(current_nodes(1,3),1));
    
    centroid_y = (1/3)*(node_list(current_nodes(1,1),2) + ...
        node_list(current_nodes(1,2),2) + node_list(current_nodes(1,3),2));
        
    % CONSTRUCT K CONTRIBUTION FOR CURRENT TRIANGLE WITHOUT BOUNDARY 
    % CONDITIONS

    for r = 1:3
        for s = r:3
            % Determine if nodes r and s are on boundary, neglecting check 
            % for exit
            node_r = current_nodes(1,r);
            node_s = current_nodes(1,s);
            is_r_on_sides = sum(wall_edge_nodes == node_r);
            is_r_on_exit = sum(exit_edge_nodes == node_r);
            is_s_on_sides = sum(wall_edge_nodes == node_s);
            is_s_on_exit = sum(exit_edge_nodes == node_s);

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
            
            linear = triangle_area*(-k0^2)*trial_r*trial_s;
            
            LHS = laplacian + linear;

            if is_r_on_sides == 0 && is_s_on_sides == 0

                 K(node_r,node_s) = K(node_r,node_s) + LHS;

                if r ~= s
                    K(node_s,node_r) = K(node_s,node_r) + LHS;   
                end

            elseif is_r_on_sides == 1 && is_s_on_sides == 0

                K(node_s,node_r) = K(node_s,node_r) + LHS;              

            elseif is_r_on_sides == 0 && is_s_on_sides == 1

                K(node_r,node_s) = K(node_r,node_s) + LHS;

            end
        end
    end

    % CONSTRUCT F CONTRIBUTION WITH PORT AND ABSORBING BC

    for r = 1:3
        node_rF = current_nodes(1,r);
        
        is_rF_on_boundary = sum(wall_edge_nodes == node_rF) + ...
            sum(exit_edge_nodes == node_rF) + ...
            sum(port_edge_nodes == node_rF);
        is_rF_on_exit_boundary = sum(exit_edge_nodes == node_rF);
        is_rF_on_port_boundary = sum(port_edge_nodes == node_rF);

        trial_rF = linear_basis_coefficients(1,r)+...
                linear_basis_coefficients(2,r)*centroid_x + ...
                linear_basis_coefficients(3,r)*centroid_y;
            
        % Estimate integral using one point gaussian quadrature 
        increment = triangle_area*current_term*trial_rF;
        
        if is_rF_on_boundary == 1
            
            % THESE INTEGRALS ARE INCORRECT
            if is_rF_on_exit_boundary == 1 % ABSORBING BOUNDARY CONDITION 

                F(node_rF,1) = F(node_rF,1) + increment;

                % First-order
                %K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*trial_rF; 
                % Second-order
                K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*(1-0.5*(pi*m/width)^2/k0^2)*trial_rF;
            end

            if is_rF_on_port_boundary == 1 % PORT BOUNDARY CONDITION 



                F(node_rF,1) = F(node_rF,1) + increment - triangle_area*2*1i*k0*init_E*sin(pi*m*centroid_y/width)*exp(1i*k0*centroid_x);

                K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*trial_rF;
            end
        
        elseif is_rF_on_boundary == 0
            F(node_rF,1) = F(node_rF,1) + increment;
        end
    end
    
    % Add side PEC boundary conditions

    for j = 1:length(wall_edge_nodes)
        BC_current_node = wall_edge_nodes(j);
        K(BC_current_node,BC_current_node) = 1;
        F(BC_current_node,1) = 0;
    end

    % Add side periodic boundary conditions (IN PROGRESS - MIGHT NEED
    % CONTINUITY IN VALUE, BELOW, and CONTINUITY IN THE FIRST DERIVATIVE IN
    % MAIN LOOP)

%     for j = 1:length(bottom_edge_nodes)
%         BC_current_node_bottom = bottom_edge_nodes(j);
%         BC_current_node_top = top_edge_nodes(j);
%         % Continuity in magnitude
%         K(BC_current_node_bottom,BC_current_node_top) = K(BC_current_node_bottom,BC_current_node_top) - 1i*k0;        
%         K(BC_current_node_top,BC_current_node_bottom) = K(BC_current_node_top,BC_current_node_bottom) + j*k0;
%     end
    
end  

% Solve system of equations
U = K\F;

% Initial signal to use for analysis
E_initial = init_E.*sin(pi*m*node_list(:,2)/width).*exp(1i*k0*node_list(:,1));

% Calculate phase of wave to test for propagation (sawtooth = good!)
phase = atan2(imag(U),real(U));

% Calculate code solution "in time" and analytic time solution
rF_cycles = 20;
[E_time,E_analytic_time] = timeAnalysis(U,node_list,omega,beta,width,init_E,rF_cycles);

num_pts = 200;
num_sweeps = 1;
windowed = 0;
% FFT analysis of solution
[k_axis,spectrum,k_parallel,R,k_calc,power] = FFTanalysis(k0,node_list,U,len,width,num_pts,'natural',windowed);
figure
plot(k_axis,spectrum)
title('FFT of E_z')
xlabel('k (m^{-1})')
ylabel('|fft(E_z)|')
% FFT analysis of initial signal
[k_vec_init,k_spectrum_init] = FFTanalysis(k0,node_list,E_initial,len,width,num_pts,'natural',windowed);
figure
plot(k_vec_init,k_spectrum_init)
title('FFT of E_z^{inc}')
xlabel('k (m^{-1})')
ylabel('|fft(E_z^{inc})|')

%----------------------------
%   PLOTTING
%----------------------------

% solution magnitude
if mag_E == 1
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),abs(U),...
        'edgecolor','k','facecolor','interp');
    view(2),axis image,colorbar
    title('Magnitude of E_z')
    xlabel('X')
    ylabel('Y')
end

% solution real component
if real_E == 1
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(U),...
        'edgecolor','k','facecolor','interp');
    view(2),axis image,colorbar
    title('Real part of E_z')
    xlabel('X')
    ylabel('Y')
end

% solution imaginary component
if imag_E == 1
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),imag(U),...
        'edgecolor','k','facecolor','interp');
    view(2),axis image,colorbar
    title('Imaginary part of E_z')
    xlabel('X')
    ylabel('Y')
end

% Plot phase calculation for solution
if phase_E == 1
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),phase,...
        'edgecolor','k','facecolor','interp');
    view(2),axis image,colorbar
    title('Phase of E_z')
    xlabel('X')
    ylabel('Y')
end

% Plot time solution movie
if time_E == 1
    figure
    peak = max(real(U));
    for j=1:length(t)
        trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(E_time(:,j)),'edgecolor','k','facecolor','interp');
        view(2),axis image, colorbar
        caxis([-peak peak])
        xlabel('X')
        ylabel('Y')
        title('Real part of E_z')
        F = getframe(gcf);
    end
end

% Plot analytic time solution movie
if analytic_time_E == 1
    figure
    peak = max(real(E_analytic_time(:,1)));
    for j=1:length(t)
        trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(E_analytic_time(:,j)),'edgecolor','k','facecolor','interp');
        view(2),axis image, colorbar
        caxis([-peak peak])
        xlabel('X')
        ylabel('Y')
        title('Analytic Solution')
        F = getframe(gcf);
    end
end

% Plot slice from a trisurf plot at the center of the waveguide (x = width/2)
if slice == 1
    num_pts = 200;
    location = width/2;
    method = 'natural';
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),U,xy(:,1),xy(:,2),method);
    figure
    plot(xy(:,1),z);
end

% Plot MOOSE solution and calculate RMS error
if MOOSE_comparison == 1
    sortMOOSEoutput
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_MOOSE,'edgecolor','k','facecolor','interp')
    view(2),axis image,colorbar
    title('MOOSE solution');
    xlabel('X')
    ylabel('Y')
    
    RMS_MOOSE_error = sqrt(sum((U_MOOSE - U).^2)/length(U))
end
        