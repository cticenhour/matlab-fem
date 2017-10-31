% buildKgraddiv.m
% Builds K components for weak form of -grad(div U) for nodes that are free 
% REQUIRES:
% K = system (KU=F) matrix
% triangle_list = element connectivity information
% node_list = node coordinates
% nonfree_edge_nodes = nodes which are not free (defined directly by some
%                       BC)
% k_x = wavenumber of direction into the page for 2D, three component 
%       problem (assumed to be x)
% component = component being modeled (0,1,2 translates to X,Y,Z respectively)
% OUTPUTS:
% K = updated system matrix

function K = buildKgraddiv(K,triangle_list,node_list,defined_edge_nodes,component,k_x)

    shift = size(node_list,1);
    
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

                % Determine gradients in YZ from coefficients of basis fxn            
                [trial_r,gradZ_r,gradY_r] = basis(current_coords,r,centroid);
                [trial_s,~,~] = basis(current_coords,s,centroid);

                % Estimate term using one point gaussian quadrature for X,
                % Y, or Z equation
                if component == 0
                    term_x = -triangle_area*(k_x^2)*trial_r*trial_s;
                    term_y = -triangle_area*k_x*gradY_r*trial_s;
                    term_z = -triangle_area*k_x*gradZ_r*trial_s;
                elseif component == 1
                    term_x = -triangle_area*k_x*gradY_r*trial_s;
                    term_y = triangle_area*gradY_r*gradY_s;
                    term_z = triangle_area*gradZ_r*gradY_s;
                elseif component == 2
                    term_x = -triangle_area*k_x*gradZ_r*trial_s;
                    term_y = triangle_area*gradY_r*gradZ_s;
                    term_z = triangle_area*gradZ_r*gradZ_s;
                else                 
                    errmsg('Set an appropriate component number! (0,1,2)');
                end

                if is_r_defined == 0 && is_s_defined == 0

                     K(node_r,node_s) = K(node_r,node_s) + term_x;
                     K(node_r+shift,node_s+shift) = K(node_r+shift,node_s+shift) + term_y;
                     K(node_r+2*shift,node_s+2*shift) = K(node_r+2*shift,node_s+2*shift) + term_z;

                    if r ~= s
                        K(node_s,node_r) = K(node_s,node_r) + term_x;
                        K(node_s+shift,node_r+shift) = K(node_s+shift,node_r+shift) + term_y;
                        K(node_s+2*shift,node_r+2*shift) = K(node_s+2*shift,node_r+2*shift) + term_z;
                    end

                elseif is_r_defined == 1 && is_s_defined == 0

                    K(node_s,node_r) = K(node_s,node_r) + term_x;
                    K(node_s+shift,node_r+shift) = K(node_s+shift,node_r+shift) + term_y;
                    K(node_s+2*shift,node_r+2*shift) = K(node_s+2*shift,node_r+2*shift) + term_z;              

                elseif is_r_defined == 0 && is_s_defined == 1

                    K(node_r,node_s) = K(node_r,node_s) + term_x;
                    K(node_r+shift,node_s+shift) = K(node_r+shift,node_s+shift) + term_y;
                    K(node_r+2*shift,node_s+2*shift) = K(node_r+2*shift,node_s+2*shift) + term_z;

                end
            end
        end
    end
end