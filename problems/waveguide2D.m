% Finite element solution of the scalar wave equation
% in a waveguide with perfectly conducting walls
% 
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% December 2016
% last update: March 10, 2017

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
    
% Initialize parts of system KU=F, where U is solution vector
K = zeros(num_nodes,num_nodes);
F = zeros(num_nodes,1);

% Weak Form 
% Term 1 - Laplacian (\nabla^2 \dot Field)
K = buildKlaplacian(K,triangle_list,node_list,wall_edge_nodes,-1);
% Term 2 - Coefficient * Field
K = buildKcoeff(K,triangle_list,node_list,wall_edge_nodes,k0);
% Right hand source term
F = buildFsource(F,triangle_list,node_list,wall_edge_nodes,source);

% Boundary Conditions
% Port BC at z = 0
[K,F] = portBC(K,F,triangle_list,node_list,port_edges, port_edge_nodes,k0);
% Absorbing BC at z = 80
[K,F] = absorbingBC(K,F,triangle_list,node_list,exit_edges,exit_edge_nodes,k0);
% PEC (E = 0) condition on walls
[K,F] = pecBC(K,F,wall_edge_nodes);

% Solve system of equations
U = K\F;

%============================
%       ANALYSIS
%============================

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

%============================
%       PLOTTING
%============================

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
    plot(xy(:,1),z,'-*');
    xlabel('Z')
    ylabel('Real(E)')
    
    analytic_real = init_E*sin(pi/width*(width/2))*cos(-sqrt(k0^2 - (pi/width)^2)*xy(:,1));
    hold on
    plot(xy(:,1),analytic_real,'-*')
    hold off
    legend('MATLAB','analytic')
    
%     abs_err_real = abs(z - analytic_real);
%     plot(xy(:,1),abs_err_real,'-*')
%     xlabel('Z')
%     ylabel('Abs. Error, Real')
    
end

% Plot slice from an imag trisurf plot at the center of the waveguide (y = width/2)
if slice_imag == 1
    num_pts = 200;
    location = width/2;
    method = 'natural';
    xy = [linspace(0,len,num_pts)', ones(num_pts,1).*(location)];
    z = griddata(node_list(:,1),node_list(:,2),imag(U),xy(:,1),xy(:,2),method);
    figure
    plot(xy(:,1),z,'-*');
    xlabel('Z')
    ylabel('Imag(E)')
    
    analytic_imag = init_E*sin(pi/width*(width/2))*sin(-sqrt(k0^2 - (pi/width)^2)*xy(:,1));
    hold on
    plot(xy(:,1),analytic_imag,'-*')
    hold off
    legend('MATLAB','analytic')

%     abs_err_imag = abs(z - analytic_imag);
%     plot(xy(:,1),abs_err_imag,'-*')
%     xlabel('Z')
%     ylabel('Abs. Error, Imaginary')
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
        