function v_sel = get_var_sel(self, variables)
%GET_VAR_SEL Get a selector for the specified variables

v_sel = false(self.mesh.n_node, length(self.sys_variables));
v_sel(:, ismember(self.sys_variables, variables)) = true;
v_sel = reshape(v_sel, [], 1);
end

