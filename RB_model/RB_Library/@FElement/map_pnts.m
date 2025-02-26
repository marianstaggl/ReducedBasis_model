function new_pos = map_pnts( self, new_support, old_pos )
%MAP_PNT this function maps the nodes of the reference element
%   the support nodes of the reference element are mapped onto the real
%   element and all the inner nodes are shifted accordingly. it behaves
%   like the function get shifted nodes but takes as old corner always the
%   ones from the reference element.

% sort the corner nodes of the linear base element
if size(new_support,1) == 4
    new_support = self.get_sorted_corners(new_support);
end

% get the support nodes of the reference element
new_pos = self.get_b(old_pos)'*new_support;

end

