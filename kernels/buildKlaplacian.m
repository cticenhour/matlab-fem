% buildKlaplacian.m
% Builds K components for weak form of laplacian with coefficient for nodes that are free
%  Regular form: Coeff * (-U_xx - U_yy)
% REQUIRES:
% K = system (KU=F) matrix
% triangle_list = element connectivity information
% node_list = node coordinates
% defined_edge_nodes = nodes which are not free (defined directly by some
%                       BC)
% coeff = coefficient for laplacian
% OUTPUTS:
% K = updated system matrix

function K = buildKlaplacian(K,triangle_list,node_list,defined_edge_nodes,coeff)

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
                [~,gradX_r,gradY_r] = basis(current_coords,r,centroid);
                [~,gradX_s,gradY_s] = basis(current_coords,s,centroid);

                % Estimate weak form term using one point gaussian quadrature
                term = coeff*triangle_area*(gradX_r*gradX_s + gradY_r*gradY_s);

                if is_r_defined == 0 && is_s_defined == 0

                     K(node_r,node_s) = K(node_r,node_s) + term;

                    if r ~= s
                        K(node_s,node_r) = K(node_s,node_r) + term;   
                    end

                elseif is_r_defined == 1 && is_s_defined == 0

                    K(node_s,node_r) = K(node_s,node_r) + term;              

                elseif is_r_defined == 0 && is_s_defined == 1

                    K(node_r,node_s) = K(node_r,node_s) + term;

                end
            end
        end
    end
end