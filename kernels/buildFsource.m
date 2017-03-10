% buildFsource.m
% Builds F components triangle-by-triangle for interior nodes only
% REQUIRES:
% F = system (KU=F) right hand side vector
% triangle_list = element connectivity information
% node_list = node coordinates
% edge_nodes = list of nodes on the edge of the domain
% source = value of source term
% OUTPUTS:
% F = updated system right hand side vector

function F = buildFsource(F,triangle_list,node_list,edge_nodes,source)

    for i = 1:size(triangle_list,1)

        current_nodes = triangle_list(i,:);
        current_coords = [node_list(current_nodes(1,1),:);...
                          node_list(current_nodes(1,2),:);...
                          node_list(current_nodes(1,3),:);];

        [triangle_area,centroid] = elementdata(current_coords);

        for r = 1:3
            node_rF = current_nodes(1,r);

            is_rF_on_boundary = sum(edge_nodes == node_rF);

            trial_rF = basis(current_coords,r,centroid);

            % Estimate integral using one point gaussian quadrature 
            increment = triangle_area*source*trial_rF;

            if is_rF_on_boundary == 0
                F(node_rF,1) = F(node_rF,1) + increment;
            end
        end
    end
end