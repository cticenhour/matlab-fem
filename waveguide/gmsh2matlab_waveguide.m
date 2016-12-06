% gmsh2matlab.m
% Script to convert Gmsh nodes, elements, and edges (as taken from MSH
% file) to usable Matlab arrays and vectors
%
% PRE-REQS: Create three plaintext files from a given MSH file:
%       1) nodes.txt -- all node data points, from number-of-nodes line to
%                       $EndNodes
%       2) elements.txt -- all 2D element connectivity information, from
%                          beginning of elm-type(second entry in each 
%                          element line) = 2 to $EndElements
%       3) edges.txt -- all 1D element connectivity information, from
%                       number-of-elements line to beginning of elm-type=2
%       *) See http://gmsh.info/doc/texinfo/gmsh.html for MSH ASCII file
%          format information, if required.
%
% NOTE: Currently configured for one boundary edge (i.e. a uniform BC in
% an FEM code)
%
% Casey Icenhour
% November 29, 2016

% clear all
% close all

load nodes.txt
load elements.txt
load sides.txt
load top.txt
load bottom.txt

num_triangles = length(elements);
num_nodes = length(nodes);
num_side_edges = size(sides,1);
num_top_edges = size(top,1);
num_bottom_edges = size(bottom,1);

triangle_list = zeros(num_triangles,3);
node_list = zeros(num_nodes,2);
side_edges = zeros(num_side_edges,2);
top_edges = zeros(num_top_edges,2);
bottom_edges = zeros(num_bottom_edges,2);

for i=1:num_nodes
    node_list(i,1:2) = [nodes(i,2) nodes(i,3)];
end

for i=1:num_triangles
    triangle_list(i,1:3) = [elements(i,6) elements(i,7) elements(i,8)];
end

for i=1:num_side_edges
    side_edges(i,1:2) = [sides(i,6) sides(i,7)];
end

for i=1:num_top_edges
    top_edges(i,1:2) = [top(i,6) top(i,7)];
end

for i=1:num_bottom_edges
    bottom_edges(i,1:2) = [bottom(i,6) bottom(i,7)];
end

side_edge_nodes = unique(side_edges);
top_edge_nodes = unique(top_edges);
bottom_edge_nodes = unique(bottom_edges);

