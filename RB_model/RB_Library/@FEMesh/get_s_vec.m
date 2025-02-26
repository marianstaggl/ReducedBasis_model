function svector = get_s_vec(line)
%SET_S_VEC calculate the s vector of the line
%   the s vector is calculated from the node positions and returns a length
%   vector starting from 0 up till the final length of the given line

% calc ds vector
ds = sqrt(diff(line.pos(:,1)).^2 + diff(line.pos(:,2)).^2);

% set the s vector and the length
svector = [0; cumsum(ds)];
end

