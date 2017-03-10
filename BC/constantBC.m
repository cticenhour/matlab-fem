% constantBC.m
% Add constant boundary condition to surfaces
% REQUIRES:
% K = system matrix
% F = system right hand side vector
% pec_edge_nodes = edge nodes where BC should be defined
% constant = constant for BC
% OUTPUTS:
% K = updated system matrix
% F = updated right hand side vector

function [K,F] = constantBC(K,F,pec_edge_nodes,constant)

    for j = 1:length(pec_edge_nodes)
        BC_current_node = pec_edge_nodes(j);
        K(BC_current_node,BC_current_node) = 1;
        F(BC_current_node,1) = constant;
    end

end