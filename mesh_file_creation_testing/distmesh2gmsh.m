% distmesh2gmsh.m
% DESCRIPTION: MATLAB script to create Gmsh ASCII format output files from 
% distMesh mesh generator output data. Based on distMesh 1.1 and 
% Gmsh 2.14.1. See NOTES below for capability specifics and usage notes.

% Created by: Casey Icenhour
% 11/22/2016
% See https://github.com/cticenhour/matlab-fem for licensing info

% distMesh belongs to Per-Olof Persson and Gilbert Strang (see 
%       http://persson.berkeley.edu/distmesh/ for more info)
% GMSH belongs to Christophe Geuzaine and Jean-François Remacle
%       (see http://gmsh.info/ for more info)

%---------------------------------

% Query output file name from user
name_query = 'What is the desired output file name (without the ".msh")? ';
fileName = input(name_query,'s');
fileName = strcat(fileName,'.msh');

% Create/open file and add header information
fileID = fopen(fileName,'wt');
fprintf(fileID,'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n');

% Create physical names (only boundary ID for now)
boundary_ID = '"edge"\n';
boundary_dim = 1;
boundary_tag = 1;
num_phys_names = 1;
fprintf(fileID,'%d\n%d %d ',num_phys_names,boundary_dim,boundary_tag);
fprintf(fileID,boundary_ID);
fprintf(fileID,'$EndPhysicalNames\n$Nodes\n');

% Write node positions
num_nodes = length(node_list);
fprintf(fileID,'%d\n',num_nodes);
for i=1:num_nodes
   x = node_list(i,1);
   y = node_list(i,2);
   z = 0;
   if abs(x) < 1e-10  % "Negative zero" correction for MATLAB XY data
       x = 0;
   end
   if abs(y) < 1e-10
       y = 0;
   end
   fprintf(fileID,'%d %1.16f %1.16f %d\n', i, x, y, z);
end
fprintf(fileID,'$EndNodes\n$Elements\n');

% Write element and connectivity data (with only 1 boundary ID)
num_elements = length(triangle_list);
fprintf(fileID,'%d\n',num_elements);
for i=1:num_elements
    % SEE GMSH DOCUMENTATION FOR THIS PART
end
fprintf(fileID,'$EndElements');

% Close and write file
fclose(fileID);
