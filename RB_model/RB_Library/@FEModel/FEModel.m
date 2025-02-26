classdef FEModel < handle
    %FEMODEL This class contains everything to solve a FE-System
    %   This class contains the mesh with its reference and parametric
    %   elements and functions to assemble the global stiffness matrix of a
    %   given system of equations. the elements are defined as shown in the
    %   sketch below eventhough the can vary if you choose a higher element
    %   order
    %
    %   3_ _ _ _ _ _4
    %   |           |
    %   |  g1   g2  |
    %   |           |
    %   |  g4   g3  |
    %   |_ _ _ _ _ _|
    %   1           2
    
    properties
        i_DoF
        sys_equations = '';
        sys_variables
        mesh = FEMesh.empty(); % mesh data
        boundary; wall_dist; wall_ids;% boundary data and wall distance
    end
    
    properties (Access = public)
        % neumann/dirichlet boundary selector and inner degrees of freedom
        dir_bounds, neu_bounds, rob_bounds
        
        % savs mapped derivatives to speed up calculations
        weights, b, dbdx, dbdy, d2bdx2, d2bdy2
        inv_jac, dbdx_ref, dbdy_ref
    end
    
    methods

        function self = FEModel(sys_equations)
            %FEMODEL Construct an instance of this class
            %   the topology of the grid is given by the elements and their
            %   connections and the boundary conditions are also passed
            
            % set the systems equations
            self.sys_equations = sys_equations;
        end
    end
    
    methods (Static)
        
        % get the default options for the solver
        opts = get_opts(varargin)
    end
end

