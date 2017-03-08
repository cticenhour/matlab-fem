% Sample removed from waveguide2D.m file for backup purposes

% Generate 2D mesh using DistMesh (see 
% http://persson.berkeley.edu/distmesh for more info)
% KEY NOTES FOR DISTMESH
% node_list = node positions -- N-by-2 array contains x,y coordinates
%                               for each of the N nodes
% triangle_list = triangle indices -- row associated with each triangle 
%                                     has 3 integer entries to specify 
%                                     node numbers corresponding to 
%                                     rows in node_list

len = 80;
width = 10;

initial_edge_width = min(len,width)/5;

geo_dist_func = @(p) drectangle(p,0,len,0,width);

bounds = [0,0;len,width];
important_pts = [0,0;len,0;0,width;len,width];

% Note: @huniform refers to uniform mesh
% See http://persson.berkeley.edu/distmesh/ for usage info
[node_list, triangle_list] = distmesh2d(geo_dist_func,@huniform,...
    initial_edge_width,bounds,important_pts);

% close distmesh rendering of mesh
close(1)

boundary_edges = boundedges(node_list,triangle_list);
edge_nodes = unique(boundary_edges);

edge_coords = node_list(edge_nodes(:),:);

port_logic = ismember(edge_coords(:,1),0);
exit_logic = ismember(edge_coords(:,1),len);
top_side_logic = ismember(edge_coords(:,2),width);
bottom_side_logic = ismember(edge_coords(:,2),0);

exit_edge_nodes = nonzeros(exit_logic.*edge_nodes);
port_edge_nodes = nonzeros(port_logic.*edge_nodes);
top_edge_nodes = nonzeros(top_side_logic.*edge_nodes);
bottom_edge_nodes = nonzeros(bottom_side_logic.*edge_nodes);

% set duplicate side corners to null values, removing them from top and
% bottom
top_edge_nodes(1) = [];
top_edge_nodes(length(top_edge_nodes)) = [];

bottom_edge_nodes(1) = [];
bottom_edge_nodes(length(bottom_edge_nodes)) = [];

wall_edge_nodes = [top_edge_nodes; bottom_edge_nodes];


num_nodes = size(node_list,1);
num_triangles = size(triangle_list,1);