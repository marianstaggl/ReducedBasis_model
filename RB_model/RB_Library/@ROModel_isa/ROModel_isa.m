classdef ROModel_isa < handle
    %ROMODEL This class sets up a reduced order model
    %   This class defines a reduced order model for the isentropic navier
    %   stokes equations. To initialize the model, an FEModel a selector
    %   for the dirichlet boudnaries and the options are passed
    
    properties
        root_folder
        visc_variables
        fe_model, def_geo % default node positions, dbdx, dbdy, and weights
        bound_nodes % selector for dirichlet bounds
        diff_sys, conv_sys
        
        base_name
        prim_shape, prim_test, visc_shape
        prim_bound, visc_bound
    end

    properties (Dependent)
        n_bound
        n_shape
        n_test
    end
    
    methods
        function self = ROModel_isa(root_folder)
            %ROMODEL Construct an instance of this class
            %   as arguments, a root folder and the dirichlet boundaries
            %   are given
            self.root_folder = root_folder;
            self.visc_variables = {'visc'};
            self.fe_model = self.get_fe_model(root_folder);
            self.bound_nodes = self.get_bound_nodes(root_folder);
            
            % store the default geometry within the model
            self.def_geo.def_nodes = self.fe_model.mesh.nodes;
            self.def_geo.weights = self.fe_model.weights;
            self.def_geo.dbdx = self.fe_model.dbdx;
            self.def_geo.dbdy = self.fe_model.dbdy;
        end
        
        function n_bound = get.n_bound(self)
            n_bound = sum(self.prim_bound);
        end

        function n_shape = get.n_shape(self)
            n_shape = size(self.prim_shape, 2) - self.n_bound;
        end

        function n_test = get.n_test(self)
            n_test = size(self.prim_test, 2) - self.n_bound;
        end
    end

    methods (Static)
        % get the default options for the solver
        opts = get_opts(varargin)
        
        % convert pressure field to local sound speed
        [full_sol, a_sel] = read_p2a(full_sol, opts)
        [full_sol, p_sel] = read_a2p(full_sol, opts)
    end
end

