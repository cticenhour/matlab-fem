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

clear all
close all

load nodes.txt
load elements.txt
load edges.txt

num_triangles = length(elements);
num_nodes = length(nodes);
num_edges = length(edges);

triangle_list = zeros(num_triangles,3);
node_list = zeros(num_nodes,2);
boundary_edges = zeros(num_edges,2);

for i=1:num_nodes
    node_list(i,1:2) = [nodes(i,2) nodes(i,3)];
end

for i=1:num_triangles
    triangle_list(i,1:3) = [elements(i,6) elements(i,7) elements(i,8)];
end

for i=1:num_edges
    boundary_edges(i,1:2) = [edges(i,6) edges(i,7)];
end

edge_nodes = unique(boundary_edges);

