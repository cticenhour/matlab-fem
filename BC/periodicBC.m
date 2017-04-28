% periodicBC.m
% Add a periodic boundary condition to surfaces
% REQUIRES:
% K = system matrix
% F = system right hand side vector
% node_list = list of node coordinates
% master_nodes = nodes which will be preserved 
% slave_nodes = nodes which will be replaced by translated master nodes
% translation = vector for translation of master node positions
% OUTPUTS:
% K = updated system matrix
% F = updated right hand side vector

function [K,F] = periodicBC(K,F,node_list, master_nodes,slave_nodes,translation)

    % Replace slave node positions in node_list with shifted master node
    % positions
    shifted_master_coords = node_list(master_nodes,:) + translation;
    for i = 1:length(master_nodes)
        current_coords = shifted_master_coords(i,:);
        dist = zeros(length(slave_nodes),1);
        slave_node_coords = node_list(slave_nodes,:);
        for j = 1:length(slave_node_coords(:,1))
            y1 = current_coords(1,1);
            z1 = current_coords(1,2);
            
            y2 = slave_node_coords(j,1);
            z2 = slave_node_coords(j,2);
            
            dist(j,1) = sqrt((y2 - y1)^2 + (z2 - z1)^2);
        end
        
        [~, min_index] = min(dist);
        node_list(slave_nodes(min_index),:) = current_coords;
    end
    
    % Add contributions to master nodes from surrounding elements to slave
    % nodes and vice versa, simulating a phantom row of elements outside of the simulation
    % box on either side and establishing continuity. 
    
    for i = 1:length(master_nodes)
        
        K_slave_old = K(slave_nodes(i),:);
        F_slave_old = F(slave_nodes(i),:);
        
        K(slave_nodes(i),:) = K(slave_nodes(i),:) + K(master_nodes(i),:);
        F(slave_nodes(i),:) = F(slave_nodes(i),:) + F(master_nodes(i),:);
        
         %K(master_nodes(i),:) = K(master_nodes(i),:) + K_slave_old;
         %F(master_nodes(i),:) = F(master_nodes(i),:) + F_slave_old;
    end
end