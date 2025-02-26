function plot_bound_nodes(self)
%PLOT_BOUND_NODES mark the nodes, specified as boundary

% get the total number of nodes per number of variables
n_prim = self.fe_model.mesh.n_node * length(self.fe_model.sys_variables);
n_visc = self.fe_model.mesh.n_node * length(self.visc_variables);

marked_prim = false(1, n_prim); marked_prim(self.bound_nodes.primary) = true;
marked_visc = false(1, n_visc); marked_visc(self.bound_nodes.viscos) = true;

all_marked = [marked_prim, marked_visc];
all_names = [self.fe_model.sys_variables(:)', self.visc_variables(:)'];
self.fe_model.mesh.plot_field(ones(1, n_prim + n_visc), all_names,...
    'marked_pts', all_marked);
end

