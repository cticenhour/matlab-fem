% basis.m
% Generates basis functions (for linear Lagrange basis functions) for a 
% triangular mesh in 2D and evaluates at a particular position
% REQUIRES:
% triangle_nodes = 3x2 array of node positions for triangle
% node_number = node number being the maximum-point of the basis function
%               (corresponds to row number in triangle_nodes)
% position = coordinate position for basis function evaluation
% OUTPUTS:
% result = numerical result of basis function at a particular position
% gradX = gradient of function in X direction
% gradY = gradient of function in Y direction


function [result,gradX,gradY] = basis(triangle_nodes,node_number,position)

    triangle_coord_matrix = [1, triangle_nodes(1,:);...
                             1, triangle_nodes(2,:);...
                             1, triangle_nodes(3,:);];
                         
    coefficients = inv(triangle_coord_matrix);

    const = coefficients(1,node_number);
    gradX = coefficients(2,node_number);
    gradY = coefficients(3,node_number);
    
    result = const + gradX*position(1,1) + gradY*position(1,2);


end