% elementdata.m
% Calculates area and centroid of triangular element in 2D
% REQUIRES:
% triangle_nodes = 3x2 array of node positions for triangle
% OUTPUTS:
% area = calculated triangle area
% centroid = 1x2 array giving coordinate for centroid

function [area,centroid] = elementdata(triangle_nodes)

    triangle_coord_matrix = [1, triangle_nodes(1,:);...
                             1, triangle_nodes(2,:);...
                             1, triangle_nodes(3,:);];

    area = abs(det(triangle_coord_matrix))/2;
    
    centroid = zeros(1,2);
    
    centroid(1,1) = (1/3)*(triangle_nodes(1,1) + triangle_nodes(2,1) + ...
        triangle_nodes(3,1));
    
    centroid(1,2) = (1/3)*(triangle_nodes(1,2) + triangle_nodes(2,2) + ...
        triangle_nodes(3,2));
end
