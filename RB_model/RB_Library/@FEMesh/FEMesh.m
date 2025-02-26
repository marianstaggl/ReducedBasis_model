classdef FEMesh < handle
    %FEMESH Contains the topology of and the mapping functionalities
    %   this class contains the nodes positions and the topology saved
    %   in the elements array. the nodes, elements and edges always contain
    %   the primary values and as last column the zone id of the entities
    
    properties (Dependent)
        n_node, n_elem, n_edge, order
    end
    
    properties (SetAccess = protected)
        dim = 2; % dimension of the mesh
        par_elem = FElement.empty; % element for the calculation of the jacobian and so on
        ref_elem = FElement.empty; % element for the shape functions
        nodes = []; % array containing the nodes coordinates
        elems = []; % elements consisting of the nodes indices
        edges = []; % edges consisting of the nodes indices
        lines = struct('id', [], 'nodes', [], 'edges', [], 'elems', []); % contains the lines formed out of the edges

        % save the results of the mapping
        weights; inv_jac; % weights and inv jacobians of reference points
        dbdx_int; dbdx_ref; dbdy_int; dbdy_ref; % mappings at ref and int points
        d2bdx2_int; d2bdy2_int; % save also the mapping for higher derivatives
    end
    
    properties (Access = protected)
        % save the original linear mesh
        def_nodes % nodes of the linear base mesh
        def_elems % elements of the linear base mesh
        def_edges % edges of the linear base mesh
    end
    
    methods
        function self = FEMesh( meshfile, order )
            
            % set the order of the shape functions
            self.ref_elem = FElement(order);
            
            % set the order of the deformation of the elements
            self.par_elem = FElement(order);
            
            % read the meshfile and save into the base properties
            [self.def_nodes, self.def_edges, self.def_elems] = ...
                self.read_meshfile( meshfile ); 
            
            % set the order of the shape functions
            [self.nodes, self.edges, self.elems] = self.get_mesh( ...
                self.def_nodes, self.def_edges, self.def_elems, order );
            
            % get the mapping functions of the elements
            [self.weights, self.dbdx_int, self.dbdy_int,...
                self.d2bdx2_int, self.d2bdy2_int, self.inv_jac,...
                self.dbdx_ref, self.dbdy_ref] = self.get_map();
            
            % set the lines around the domain
            self.lines = self.get_lines( self.nodes, self.edges,...
                self.elems );
        end
        
        % define the get functions for dependent properties
        function n_node = get.n_node( self ), n_node = size(self.nodes,1); end
        
        function n_elem = get.n_elem( self ), n_elem = size(self.elems,1); end
        
        function n_edge = get.n_edge( self ), n_edge = size(self.edges,1); end
        
        function order = get.order( self ), order = self.ref_elem.order; end
    end
    
    methods (Static)
        % split the mesh up for higher order elements
        [ nodes, edges, elements ] = get_mesh( nodes, edges, elements, order )
        
        % get sorted lines from edges
        lines = get_lines( nodes, edges, elems );
        
        % get the s vector of a line
        s_vec = get_s_vec( line );
        
        % split up a line according to angle th
        lines = static_split_lines( lines, threshold )
        
        % select the lines with ginput
        ids = static_select_lines( lines, cid )
        
        % fit a circle into a line in 2D
        [xm, ym, r] = static_fit_circle( lines, cid )
        
        % create a control mesh
        c_mesh = static_get_c_mesh(b_ids, n_elem);
    end
end
