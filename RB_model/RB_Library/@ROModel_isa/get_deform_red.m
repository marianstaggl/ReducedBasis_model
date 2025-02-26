function all_red = get_deform_red(self)
%GET_DEFORM_DEIM Summary of this function goes here
%   Detailed explanation goes here

diff_elem = self.diff_sys.deim_elem;
conv_elem = self.conv_sys.deim_elem;

% get a list of all necessary elements
all_red = unique([diff_elem, conv_elem]);
end

