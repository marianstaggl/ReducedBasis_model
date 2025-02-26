function [dxdx_full, dydy_full, ww_full] = get_deform_aux(self, red_elem)
%GET_DEFORM_DXDXDYDY get the dxdx and dydy of the deformation

% get the number of interpolation points and the number of deim elems
n_elem = size(self.fe_model.weights, 1);
n_intp = size(self.fe_model.weights, 2);

% preallocate the necessary arrays
dxdx_full = zeros(n_elem, n_intp);
dydy_full = zeros(n_elem, n_intp);
ww_full = zeros(n_elem, n_intp);

% check if the nodes are dlarrays, if so, convert the rest
if isa(self.fe_model.mesh.nodes, 'dlarray')
    dxdx_full = dlarray(dxdx_full);
    dydy_full = dlarray(dydy_full); 
    ww_full = dlarray(ww_full);
end

% get the differences of the positions
ax = self.fe_model.mesh.nodes(:, 1) - self.def_geo.def_nodes(:, 1);
ay = self.fe_model.mesh.nodes(:,2) - self.def_geo.def_nodes(:,2);
ww_full(red_elem, :) = self.fe_model.weights(red_elem, :) ./ ...
    self.def_geo.weights(red_elem, :);

% calculate dxdx and dydy using a vectorized approach
elems = self.fe_model.mesh.elems(red_elem,:);
ax_full = permute(reshape(ax(elems), [], n_intp), [2, 3, 1]);
ay_full = permute(reshape(ay(elems), [], n_intp), [2, 3, 1]);
dxdx_full(red_elem, :) = 1 - permute(sum(self.fe_model.dbdx(:,:,red_elem)...
    .* ax_full, 1), [3, 2, 1]);
dydy_full(red_elem, :) = 1 - permute(sum(self.fe_model.dbdy(:,:,red_elem)...
    .* ay_full, 1), [3, 2, 1]);
end

% % loop over the elements
% for i=1:length(red_elem)
%     c_elem = red_elem(i);
%     curr_elem = self.fe_model.mesh.elems(c_elem,:); % pick the current element
%     daxdx_int = sum(self.fe_model.dbdx(:,:,c_elem).*ax(curr_elem), 1); % get the ax derivatives at the integration points
%     daydy_int = sum(self.fe_model.dbdy(:,:,c_elem).*ay(curr_elem), 1); % get the ay derivatives at the integration points
% 
%     dxdx_full(c_elem,:) = 1 - daxdx_int; % get the transformation dxdx at the integration points
%     dydy_full(c_elem,:) = 1 - daydy_int; % get the transformation dydy at the integration points
%     ww_full(c_elem,:) = ww_red(i, :);
% end

