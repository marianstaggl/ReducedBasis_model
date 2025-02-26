function [n1_vecs, n2_vecs] = get_line_n(self, line_no)
%GET_LINE_NORMALS Get the line normals for boundary treatment.
%   the n1_vecs are normal to the wall, the n2_vecs are parallel to the
%   wall. both are used to get d(vel)/d(n1_vec) the wall shear stresses at
%   the wall points

% pick the current line to get the node ids
c_line = self.lines(line_no);

% calculate the line normals (at the center of the pieces)
n2_vecs = normr(c_line.pos(2:end,:) - c_line.pos(1:end-1,:)); % p2p vecs
n1_vecs = normr([-n2_vecs(:,2), n2_vecs(:,1)]); % normalized normals

% flip the normals if necessary to point them inwards
end

