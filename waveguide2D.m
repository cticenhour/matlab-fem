% Finite element solution of the scalar wave equation
% in a waveguide with perfectly conducting walls
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% December 2016
% last update: March 8, 2017

%=============================
% SWITCHES AND PLOTTING OPTIONS
%=============================
    
filename = 'waveguide.msh';

MOOSE_comparison = 0;   % requires output data CSV file from MOOSE

% Plotting switches
surface_plots = 1;
phase_E = 0;
time_E = 0;
analytic_time_E = 0;
slice_real = 0;
slice_imag = 0;
fft = 0;

%=============================
% IMPORTANT CONSTANTS
%=============================

width = 10;
len = 80;

m = 1;

init_E = 1;

omega = 2*pi*20e6;

mu0 = 4*pi*1e-7;

eps0 = 8.854e-12; % F/m

c = 3e8;

Z0 = mu0*c;

k0 = omega/c;

k = k0*(1+1i*0);

beta = sqrt(k^2 - (pi/width)^2);

source = 0;

%=============================

% Create usable connectivity info from GMSH mesh file

[node_list,triangle_list,boundary_edges,boundary_names] = ...
                                    gmsh2matlab2d(filename);
% Create edge nodes and edges arrays, labeled with appropriate names
total_bounds = length(boundary_names);
for i = 1:total_bounds
    eval([boundary_names{i},'_edge_nodes = nonzeros(unique(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,1:2)));'])
    eval([boundary_names{i},'_edges = [nonzeros(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,1)), nonzeros(ismember(boundary_edges(:,3),',num2str(i),').*boundary_edges(:,2))];'])
end

% Clean up task - Remove duplicate corners on exit and port from walls
top_edge_nodes(ismember(top_edge_nodes,intersect(top_edge_nodes,port_edge_nodes))) = [];
top_edge_nodes(ismember(top_edge_nodes,intersect(top_edge_nodes,exit_edge_nodes))) = [];

bottom_edge_nodes(ismember(bottom_edge_nodes,intersect(bottom_edge_nodes,port_edge_nodes))) = [];
bottom_edge_nodes(ismember(bottom_edge_nodes,intersect(bottom_edge_nodes,exit_edge_nodes))) = [];

% Create other needed things
wall_edge_nodes = [top_edge_nodes; bottom_edge_nodes];
num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);
    
% Initialize parts of system KU=F, where u is solution vector
K = zeros(num_nodes,num_nodes);
F = zeros(num_nodes,1);

% Applying element-oriented FEM algorithm to construct K and F
for i = 1:num_triangles
    
    current_nodes = triangle_list(i,:);
    current_coords = [node_list(current_nodes(1,1),:);...
                      node_list(current_nodes(1,2),:);...
                      node_list(current_nodes(1,3),:);];
    
    [triangle_area,centroid] = elementdata(current_coords);
    centroid_x = centroid(1,1);
    centroid_y = centroid(1,2);
        
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
            [trial_r,gradX_r,gradY_r] = basis(current_coords,r,centroid);
            [trial_s,gradX_s,gradY_s] = basis(current_coords,s,centroid);
            
            % Estimate first term using one point gaussian quadrature
            laplacian = -triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);
            
            % Estimate second term using one point gaussian quadrature
            linear = triangle_area*(k0^2)*trial_r*trial_s;
            
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
        
%         is_rF_on_boundary = sum(wall_edge_nodes == node_rF) + ...
%             sum(exit_edge_nodes == node_rF) + ...
%             sum(port_edge_nodes == node_rF);
        is_rF_on_exit_boundary = sum(exit_edge_nodes == node_rF);
        is_rF_on_port_boundary = sum(port_edge_nodes == node_rF);

        trial_rF = basis(current_coords,r,centroid);
            
        % Estimate integral using one point gaussian quadrature 
        increment = triangle_area*source*trial_rF;
        
        if is_rF_on_exit_boundary == 1 || is_rF_on_port_boundary == 1
            
            if is_rF_on_exit_boundary == 1 % ABSORBING BOUNDARY CONDITION 
                
                F(node_rF,1) = F(node_rF,1) + increment;             

                % First-order - OLD INCORRECT INTEGRALS
                %K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*trial_rF; 
                % Second-order - OLD
                %K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*(1-0.5*(pi*m/width)^2/k0^2)*trial_rF;
                
                % INTEGRATIONS ARE PERFORMED USING THE COMPOSITE SIMPSON'S
                % RULE
                
                % Get other endpoint of boundary edge for node_rF              
                adj_nodes = unique(nonzeros(exit_edges.*sum(ismember(exit_edges,node_rF),2)));
                node2 = adj_nodes(ismember(adj_nodes,current_nodes) &~ ismember(adj_nodes,node_rF));
                
                if sum(node2) ~= 0 % If node_rF is on a triangle with a full edge on the boundary
                    halfway_y = 0.5*(node_list(node2,2) + node_list(node_rF,2));

                    % Calculate values for integration based on relative
                    % node positions
                    if node_list(node_rF,2) > node_list(node2,2)
                        trial_a = basis(current_coords,r,node_list(node2,:));
                        trial_half = basis(current_coords,r,[node_list(node2,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node_rF,:)); 
                            
                        b_minus_a = node_list(node_rF,2) - node_list(node2,2);
                    else
                        
                        trial_a = basis(current_coords,r,node_list(node_rF,:));
                        trial_half = basis(current_coords,r,[node_list(node_rF,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node2,:));   
                            
                        b_minus_a = node_list(node2,2) - node_list(node_rF,2);
                    end
                else % If node_rF is simply a point on the edge
                    trial_a = 0;        % Sets integration to zero - 
                    trial_half = 0;     % cannot integrate over a point and 
                    trial_b = 0;        % integration is covered in
                    b_minus_a = 0;      % triangles with edges on the
                                        % boundary

                end
                
                K(node_rF,node_rF) = K(node_rF,node_rF) - (b_minus_a/6)*1i*sqrt(k0^2 - (pi*m/width)^2)*(trial_a + 4*trial_half + trial_b);
                %K(node_rF,node_rF) = K(node_rF,node_rF) - (b_minus_a/6)*1i*sqrt(k0^2 - (pi*m/width)^2)*(1-0.5*(pi*m/width)^2/sqrt(k0^2 - (pi*m/width)^2)^2)*(trial_a + 4*trial_half + trial_b);
            end

            if is_rF_on_port_boundary == 1 % PORT BOUNDARY CONDITION 

                % INTEGRATIONS ARE PERFORMED USING THE COMPOSITE SIMPSON'S
                % RULE
                
                % Get other endpoint of boundary edge for node_rF              
                adj_nodes = unique(nonzeros(port_edges.*sum(ismember(port_edges,node_rF),2)));
                node2 = adj_nodes(ismember(adj_nodes,current_nodes) &~ ismember(adj_nodes,node_rF));
                
                if sum(node2) ~= 0 % If node_rF is on a triangle with a full edge on the boundary
                    halfway_y = 0.5*(node_list(node2,2) + node_list(node_rF,2));

                    % Calculate values for integration based on relative
                    % node positions
                    if node_list(node_rF,2) > node_list(node2,2)
                        trial_a = basis(current_coords,r,node_list(node2,:));
                        trial_half = basis(current_coords,r,[node_list(node2,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node_rF,:));                        

                        inc_a = sin(pi*m*node_list(node2,2)/width)*exp(1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node2,1));
                        inc_half = sin(pi*m*halfway_y/width)*exp(1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node2,1));
                        inc_b = sin(pi*m*node_list(node_rF,2)/width)*exp(1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node_rF,1));
                        
                        b_minus_a = node_list(node_rF,2) - node_list(node2,2);
                    else
                        trial_a = basis(current_coords,r,node_list(node_rF,:));
                        trial_half = basis(current_coords,r,[node_list(node_rF,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node2,:));

                        inc_a = sin(pi*m*node_list(node_rF,2)/width)*exp(-1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node_rF,1));
                        inc_half = sin(pi*m*halfway_y/width)*exp(-1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node_rF,1));
                        inc_b = sin(pi*m*node_list(node2,2)/width)*exp(-1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(node2,1));
                        
                        b_minus_a = node_list(node2,2) - node_list(node_rF,2);
                    end
                    
                else % If node_rF is simply a point on the edge
                    trial_a = 0;        % Sets integration to zero - 
                    trial_half = 0;     % cannot integrate over a point and 
                    trial_b = 0;        % integration is covered in
                    inc_a = 0;          % triangles with edges on the
                    inc_half = 0;       % boundary
                    inc_b = 0;
                    b_minus_a = 0;
                end                    
                
                F(node_rF,1) = F(node_rF,1) + increment - (b_minus_a/6)*1i*sqrt(k0^2 - (pi*m/width)^2)*init_E*(inc_a + 4*inc_half + inc_b);

                K(node_rF,node_rF) = K(node_rF,node_rF) - (b_minus_a/6)*1i*sqrt(k0^2 - (pi*m/width)^2)*(trial_a + 4*trial_half + trial_b);
                
                % OLD
%                 F(node_rF,1) = F(node_rF,1) + increment - triangle_area*2*1i*k0*init_E*sin(pi*m*centroid_y/width)*exp(1i*k0*centroid_x);
% 
%                 K(node_rF,node_rF) = K(node_rF,node_rF) - triangle_area*1i*k0*trial_rF;
            end
        
        else
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
E_initial = init_E.*sin(pi*m*node_list(:,2)/width).*exp(-1i*sqrt(k0^2 - (pi*m/width)^2)*node_list(:,1));

% Calculate phase of wave to test for propagation (sawtooth = good!)
phase = atan2(imag(U),real(U));

% Calculate code solution "in time" and analytic time solution
rF_cycles = 20;
[E_time,E_analytic_time] = timeAnalysis(U,node_list,omega,beta,width,init_E,rF_cycles);

num_pts = 200;
num_sweeps = 1;
windowed = 0;
% FFT analysis of solution
[k_axis,spectrum,k_parallel,R,k_calc,angle,R_theory_first, R_theory_second,power] = FFTanalysis(k0,node_list,U,len,width,num_pts,'natural',windowed);
% FFT analysis of initial signal
[k_vec_init,k_spectrum_init] = FFTanalysis(k0,node_list,E_initial,len,width,num_pts,'natural',windowed);

%----------------------------
%   PLOTTING
%----------------------------

if surface_plots == 1
    % solution magnitude
    figure
    subplot(3,1,1)
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),abs(U),...
        'edgecolor','k','facecolor','interp');
    view(2)
    axis image
    colorbar
    title('Magnitude of E_z')
    xlabel('Z')
    ylabel('Y')

    % solution real component
    subplot(3,1,2)
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(U),...
        'edgecolor','k','facecolor','interp');
    view(2)
    axis image
    colorbar
    title('Real part of E_z')
    xlabel('Z')
    ylabel('Y')

    % solution imaginary component
    subplot(3,1,3)
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),imag(U),...
        'edgecolor','k','facecolor','interp');
    view(2)
    axis image
    colorbar
    title('Imaginary part of E_z')
    xlabel('Z')
    ylabel('Y')
end

% Plot phase calculation for solution
if phase_E == 1
    figure
    trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),phase,...
        'edgecolor','k','facecolor','interp');
    view(2)
    axis image
    colorbar
    title('Phase of E_z')
    xlabel('Z')
    ylabel('Y')
end

% Plot time solution movie
if time_E == 1
    figure
    peak = max(real(U));
    for j=1:length(t)
        trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),real(E_time(:,j)),'edgecolor','k','facecolor','interp');
        view(2)
        axis image
        colorbar
        caxis([-peak peak])
        xlabel('Z')
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
        view(2)
        axis image
        colorbar
        caxis([-peak peak])
        xlabel('Z')
        ylabel('Y')
        title('Analytic Solution')
        F = getframe(gcf);
    end
end

% Plot slice from a real trisurf plot at the center of the waveguide (y = width/2)
if slice_real == 1
    num_pts = 200;
    location = width/2;
    method = 'natural';
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),real(U),xy(:,1),xy(:,2),method);
    figure
%     plot(xy(:,1),z,'-*');
%     xlabel('Z')
%     ylabel('Real(E)')
    
    analytic_real = init_E*sin(pi/width*(width/2))*cos(-sqrt(k0^2 - (pi/width)^2)*xy(:,1));
%     hold on
%     plot(xy(:,1),analytic_real,'-*')
%     hold off
%     legend('MATLAB','analytic')
    
    abs_err_real = abs(z - analytic_real);
    plot(xy(:,1),abs_err_real,'-*')
    xlabel('Z')
    ylabel('Abs. Error, Real')
    
end

% Plot slice from an imag trisurf plot at the center of the waveguide (y = width/2)
if slice_imag == 1
    num_pts = 200;
    location = width/2;
    method = 'natural';
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),imag(U),xy(:,1),xy(:,2),method);
    figure
%     plot(xy(:,1),z,'-*');
%     xlabel('Z')
%     ylabel('Imag(E)')
    
    analytic_imag = init_E*sin(pi/width*(width/2))*sin(-sqrt(k0^2 - (pi/width)^2)*xy(:,1));
%     hold on
%     plot(xy(:,1),analytic_imag,'-*')
%     hold off
%     legend('MATLAB','analytic')

    abs_err_imag = abs(z - analytic_imag);
    plot(xy(:,1),abs_err_imag,'-*')
    xlabel('Z')
    ylabel('Abs. Error, Imaginary')
end

% FFT analysis and plotting
if fft == 1
    num_pts = 200;
    num_sweeps = 1;
    windowed = 0;
    % FFT analysis of solution
    [k_axis,spectrum,k_parallel,R,k_calc,angle,R_theory_first, R_theory_second,power] = FFTanalysis(k0,node_list,U,len,width,num_pts,'natural',windowed);
    % FFT analysis of initial signal
    [k_vec_init,k_spectrum_init] = FFTanalysis(k0,node_list,E_initial,len,width,num_pts,'natural',windowed);
    
    figure
    plot(k_axis,spectrum)
    title('FFT of E_z')
    xlabel('k (m^{-1})')
    ylabel('|fft(E_z)|')
    figure
    plot(k_vec_init,k_spectrum_init)
    title('FFT of E_z^{inc}')
    xlabel('k (m^{-1})')
    ylabel('|fft(E_z^{inc})|')
end

% Plot MOOSE solution and calculate RMS error
if MOOSE_comparison == 1
    sortMOOSEoutput
%     figure
%     trisurf(triangle_list,node_list(:,1),node_list(:,2),0*node_list(:,1),U_MOOSE_Im,'edgecolor','k','facecolor','interp')
%     view(2),axis image,colorbar
%     xlabel('X')
%     ylabel('Y')
    
    num_pts = 200;
    location = width/2;
    method = 'natural';
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),U_MOOSE_Re,xy(:,1),xy(:,2),method);
    figure(1)
    hold on
    plot(xy(:,1),z,'-*');
    hold off
    legend('MATLAB', 'analytic', 'MOOSE')
    
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),U_MOOSE_Im,xy(:,1),xy(:,2),method);
    figure(2)
    hold on
    plot(xy(:,1),z,'-*');
    hold off
    legend('MATLAB', 'analytic', 'MOOSE')
    
    %RMS_MOOSE_error = sqrt(sum((U_MOOSE - U).^2)/length(U))
end
        