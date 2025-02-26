function [deform_diff, deform_marker] = get_deform_diff(self, red_elem)
%GET_DEFORM_DIFF get the deformation function for the diffusive part

% get the derivatives of the deformation over x
[dxdx, dydy, ww] = self.get_deform_aux(red_elem);

% build the linear function for this deformation
phi_1 = reshape(dxdx .* dxdx .* ww, [], 1);
phi_2 = reshape(dydy .* dydy .* ww, [], 1);
deform_diff = reshape([phi_1, phi_2], [], 1);
deform_marker = 2; % number of marker for the deformation
end