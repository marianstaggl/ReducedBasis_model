classdef FElement < handle
    %ELEMENT This class defines a quad element
    %   the four corner nodes are defined as follows. with these nodes it
    %   is possible to split the element (to increase the order of the
    %   shapefunctions). the element object contains all functions needed
    %   to assemble the systems matrix including: the basefunction and its
    %   derivatives (always defined on the reference element).
    %
    %   3_ _ _ _ _ _4
    %   |           |
    %   |           |
    %   |           |
    %   |           |
    %   |_ _ _ _ _ _|
    %   1           2
    
    properties (SetAccess = protected)
        order % order of the elements shapefunctions
        ref_points % the points of the lagrangian support
        int_points % the integration points
        int_weights % set the integration weights
    end
    
    methods
        % constructor method for the element
        function self = FElement(order)
            
            % set the order of the element
            self.order = order;
            
            % get the reference (support) points of the lagrange polys
            self.ref_points = self.get_ref_points( order );
            
            % set the points for the gaussian quardature
            [self.int_points, self.int_weights] = self.get_int_points( order );
        end
    end
    
    methods (Static)
        
        % get a lagrange base polynomial of degree n defined by the support
        Lj = get_lagrange_base( eta, eta_support, base_idx, derivation_order );
        
        % sort the given corner nodes
        [sorted_corners, sort_idx] = get_sorted_corners( corner_nodes );
        
        % get the support nodes of the regualar element
        [node_pos, node_ids, corner_ids] = get_ref_points( order );
        
        % get the integration points of the regular element
        [points, weights] = get_int_points( order )
        
        % use a bilinear transformation to shift nodes
        nodes = get_shifted_nodes(old_supp, new_supp, nodes)
    end
end

