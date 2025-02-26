function set_deim_elem(self)
%SET_DEIM_ELEM Summary of this function goes here
%   Detailed explanation goes here

self.diff_sys.deim_elem = self.get_deform_elem(3, self.diff_sys.diff_deim);
self.conv_sys.deim_elem = self.get_deform_elem(2, self.conv_sys.conv_deim);
end

