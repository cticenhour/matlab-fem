% GMSH2MATLAB2D.m
% MATLAB script to convert Gmsh nodes, elements, and edges from output MSH
% file into usable MATLAB arrays and vectors
%
% Creator: Casey T. Icenhour
% Last Update: February 6, 2017
% 
% License: MIT License
%
% Reference: 
% 
%       C. Geuzaine and J.-F. Remacle. "Gmsh: a three-dimensional finite 
%       element mesh generator with built-in pre- and post-processing 
%       facilities." International Journal for Numerical Methods in 
%       Engineering 79(11), pp. 1309-1331, 2009.
%
% Usage:
%       gmsh2matlab2d('example.msh')
%   
% Output: 
%
%       node_list (3-column array with all node locations)
%       triangle_list (3-column array with all triangle connectivity 
%                       info - node numbers are indices in node_list)
%       boundary_edges (3-column array of boundary edge info:
%                           Column 1: first edge node index
%                           Column 2: second edge node index
%                           Column 3: boundary name ID number
%       boundary_names (cell array of physical boundary names, boundary 
%                           ID numbers refer to index of this array)

function [node_list,triangle_list,boundary_edges,boundary_names] = ...
                                                    gmsh2matlab2d(filename)

%--------------------
% Open File
%--------------------

[fileID, errmsg] = fopen(filename, 'rt');

if fileID < 0
    error(errmsg);
end

%---------------------------------------
% Get physical boundary names
%---------------------------------------

while(1) 
    line = fgetl(fileID);
    k = strfind(line,'$PhysicalNames');
    j = isempty(k);
    if j == 0
        total = str2double(fgetl(fileID));
        boundary_names = cell(total,1);
        for i=1:total
            line = fgetl(fileID);
            boundary_names{i} = line(6:length(line)-1);
        end
        break
    end
end

%---------------------------------------
% Get node locations
%---------------------------------------

while(1)
    line = fgetl(fileID);
    k = strfind(line,'$Nodes');
    j = isempty(k);
    if j == 0
       total = str2double(fgetl(fileID));
       node_list = zeros(total,3);
       for i=1:total
           line = fgetl(fileID);
           line_split = strsplit(line);
           node_list(i,1) = str2double(line_split(2));
           node_list(i,2) = str2double(line_split(3));
           node_list(i,3) = str2double(line_split(4));
       end
       break
    end
end

%---------------------------------------------
% Get edges with boundary IDs and triangles
%---------------------------------------------

while(1)
    line = fgetl(fileID);
    k = strfind(line,'$Elements');
    j = isempty(k);
    if j == 0
        total = str2double(fgetl(fileID));
        boundary_edges = zeros(1,3);
        triangle_list = zeros(1,3);
        for i = 1:total
            line = fgetl(fileID);
            line_split = strsplit(line);
            elem_type = str2double(line_split(2));
            if elem_type == 1
                boundaryID = str2double(line_split(4));
                boundary_edges(i,1) = str2double(line_split(6));
                boundary_edges(i,2) = str2double(line_split(7));
                boundary_edges(i,3) = boundaryID;
            elseif elem_type == 2
                tri_index = i - length(boundary_edges);
                triangle_list(tri_index,1) = str2double(line_split(6));
                triangle_list(tri_index,2) = str2double(line_split(7));
                triangle_list(tri_index,3) = str2double(line_split(8));
            end
        end
        break 
    end
end

fclose(fileID);

end




