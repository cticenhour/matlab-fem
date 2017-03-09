% buildKinterior.m
% Builds K components triangle-by-triangle for interior nodes only
% REQUIRES:
% K = system (KU=F) matrix
% triangle_list = element connectivity information
% node_list = node coordinates
% wall_edge_nodes = nodes located on the walls
% k = wave number for linear term
% OUTPUTS:
% K = updated system matrix

function K = buildKinterior(K,triangle_list,node_list,defined_edge_nodes,k)

    for i = 1:size(triangle_list,1)

        current_nodes = triangle_list(i,:);
        current_coords = [node_list(current_nodes(1,1),:);...
                          node_list(current_nodes(1,2),:);...
                          node_list(current_nodes(1,3),:);];

        [triangle_area,centroid] = elementdata(current_coords);

        for r = 1:3
            for s = r:3
                % Determine if nodes r and s are on boundary, neglecting check 
                % for exit
                node_r = current_nodes(1,r);
                node_s = current_nodes(1,s);
                is_r_defined = sum(defined_edge_nodes == node_r);
                is_s_defined = sum(defined_edge_nodes == node_s);

                % Determine gradients in XY from coefficients of basis fxn            
                [trial_r,gradX_r,gradY_r] = basis(current_coords,r,centroid);
                [trial_s,gradX_s,gradY_s] = basis(current_coords,s,centroid);

                % Estimate first term using one point gaussian quadrature
                laplacian = -triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);

                % Estimate second term using one point gaussian quadrature
                linear = triangle_area*(k^2)*trial_r*trial_s;

                LHS = laplacian + linear;

                if is_r_defined == 0 && is_s_defined == 0

                     K(node_r,node_s) = K(node_r,node_s) + LHS;

                    if r ~= s
                        K(node_s,node_r) = K(node_s,node_r) + LHS;   
                    end

                elseif is_r_defined == 1 && is_s_defined == 0

                    K(node_s,node_r) = K(node_s,node_r) + LHS;              

                elseif is_r_defined == 0 && is_s_defined == 1

                    K(node_r,node_s) = K(node_r,node_s) + LHS;

                end
            end
        end
    end
end