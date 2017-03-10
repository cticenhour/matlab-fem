% portBC.m
% Add Port Boundary Condition to Surfaces
% REQUIRES:
% K = system matrix
% F = right hand side vector
% triangle_list = triangle mesh connectivity info
% node_list = node coordinates list
% port_edges = edges along which portBC is defined
% port_edge_nodes = nodes where port BC is defined
% k = wave number
% OUTPUTS
% K = updated system matrix
% F = update right hand side vector

function [K,F] = portBC(K,F,triangle_list,node_list,port_edges,port_edge_nodes,k)
    
    m = 1;
    width = 10;


    for i = 1:size(triangle_list,1)

        current_nodes = triangle_list(i,:);
        current_coords = [node_list(current_nodes(1,1),:);...
                          node_list(current_nodes(1,2),:);...
                          node_list(current_nodes(1,3),:);];

        for r = 1:3
            node_rF = current_nodes(1,r);
            
            is_rF_on_port_boundary = sum(port_edge_nodes == node_rF);

            if is_rF_on_port_boundary == 1 % PORT BOUNDARY CONDITION 

                % INTEGRATIONS ARE PERFORMED USING THE COMPOSITE SIMPSON'S
                % RULE

                % Get other endpoint of boundary edge for node_rF              
                adj_nodes = unique(nonzeros(port_edges.*sum(ismember(port_edges,node_rF),2)));
                node2 = adj_nodes(ismember(adj_nodes,current_nodes) &~ ismember(adj_nodes,node_rF));

                if sum(node2) ~= 0 % If node_rF is on a triangle with a full edge on the boundary
                    halfway_y = 0.5*(node_list(node2,2) + node_list(node_rF,2));

                    % Calculate values for integration based on relative
                    % node positions
                    if node_list(node_rF,2) > node_list(node2,2)
                        trial_a = basis(current_coords,r,node_list(node2,:));
                        trial_half = basis(current_coords,r,[node_list(node2,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node_rF,:));                        

                        inc_a = sin(pi*m*node_list(node2,2)/width)*exp(1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node2,1));
                        inc_half = sin(pi*m*halfway_y/width)*exp(1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node2,1));
                        inc_b = sin(pi*m*node_list(node_rF,2)/width)*exp(1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node_rF,1));

                        b_minus_a = node_list(node_rF,2) - node_list(node2,2);
                    else
                        trial_a = basis(current_coords,r,node_list(node_rF,:));
                        trial_half = basis(current_coords,r,[node_list(node_rF,1),halfway_y]);
                        trial_b = basis(current_coords,r,node_list(node2,:));

                        inc_a = sin(pi*m*node_list(node_rF,2)/width)*exp(-1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node_rF,1));
                        inc_half = sin(pi*m*halfway_y/width)*exp(-1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node_rF,1));
                        inc_b = sin(pi*m*node_list(node2,2)/width)*exp(-1i*sqrt(k^2 - (pi*m/width)^2)*node_list(node2,1));

                        b_minus_a = node_list(node2,2) - node_list(node_rF,2);
                    end

                else % If node_rF is simply a point on the edge
                    trial_a = 0;        % Sets integration to zero - 
                    trial_half = 0;     % cannot integrate over a point and 
                    trial_b = 0;        % integration is covered in
                    inc_a = 0;          % triangles with edges on the
                    inc_half = 0;       % boundary
                    inc_b = 0;
                    b_minus_a = 0;
                end                    

                F(node_rF,1) = F(node_rF,1) - (b_minus_a/6)*1i*sqrt(k^2 - (pi*m/width)^2)*(inc_a + 4*inc_half + inc_b);

                K(node_rF,node_rF) = K(node_rF,node_rF) - (b_minus_a/6)*1i*sqrt(k^2 - (pi*m/width)^2)*(trial_a + 4*trial_half + trial_b);

            end
        end
    end
end
    