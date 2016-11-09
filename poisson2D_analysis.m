% Basic error analysis for Poisson 2D FEM code
%
% Author: Casey Icenhour
% Organization: North Carolina State University/Oak Ridge National
%                                               Laboratory
% September 2016

clear all 
close all

resolutions = 0.05:0.01:0.7;
num_elements = zeros(1,length(resolutions));
error_RMS = zeros(length(resolutions),1);


for i = 1:length(resolutions)
   [U,node_list,triangle_list] = poisson2DFXN(1,resolutions(i));
   U_analytic = analyticFXN_poisson2D(node_list,50);
   
   error_RMS(i,1) = sqrt(sum((U-U_analytic).^2)/(length(U)));
   num_elements(i,1) = length(triangle_list); 
end

close all
figure
semilogy(num_elements,error_RMS)
xlabel('Number of Elements')
ylabel('Error_{RMS}')
