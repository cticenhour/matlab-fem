% Analysis script for waveguide FEM resolution studies

clear all

factor = 0.1:0.02:1;
reflection = zeros(1,length(factor));
reflection_theory = zeros(1,length(factor));
k_calculated = zeros(1,length(factor));
elements = zeros(1,length(factor));
percentage = 0;

for analysis_index = 1:length(factor)
    clearvars -except analysis_index factor reflection elements percentage k_calculated
    command = ['gmsh -2 meshes/waveguide.geo -clscale ',num2str(factor(analysis_index))];
    s = system(command);
    if s ~= 0
        error('ERROR -gmsh command not successful!')
    end
    waveguide2D
    reflection(analysis_index) = R;
    reflection_theory(analysis_index) = R_theory_second;
    k_calculated(analysis_index) = k_calc;
    elements(analysis_index) = num_triangles;
end

figure
plot(elements,reflection.*100,'-x')
xlabel('Number of Elements')
ylabel('Reflection (%)')
title(['Min Reflection: ' num2str(min(reflection).*100) '%'])

figure
plot(elements,k_calculated./k0,'-x')
xlabel('Number of Elements')
ylabel('k_{calculated}/k_0')
