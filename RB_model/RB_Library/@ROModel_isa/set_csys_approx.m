function set_csys_approx(self)
%SET_SYS_MAT create the 3rd order tensors for the reduced system

% take the 4th order tensors and condense them into 3rd order
conv_def = self.get_deform_conv(self.conv_sys.deim_elem);
theta_conv = pinv(self.conv_sys.conv_xi) * conv_def(self.conv_sys.conv_deim); % self.conv_sys.conv_xi \ conv_def(self.conv_sys.conv_deim);
conv1_3rd = sum(self.conv_sys.tensor_1 .* reshape(theta_conv, 1, 1, 1, []), 4);
conv2_3rd = sum(self.conv_sys.tensor_2 .* reshape(theta_conv, 1, 1, 1, []), 4);

diff_def = self.get_deform_diff(self.diff_sys.deim_elem);
theta_diff = pinv(self.diff_sys.diff_xi) * diff_def(self.diff_sys.diff_deim); % self.diff_sys.diff_xi \ diff_def(self.diff_sys.diff_deim);
diff_3rd = sum(self.diff_sys.tensor .* reshape(theta_diff, 1, 1, 1 ,[]), 4);

% set the current system matrices
self.diff_sys.diff_3rd = diff_3rd;
self.conv_sys.conv1_3rd = conv1_3rd;
self.conv_sys.conv2_3rd = conv2_3rd;
end
