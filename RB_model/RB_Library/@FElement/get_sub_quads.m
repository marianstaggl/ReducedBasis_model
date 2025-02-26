function sub_quads = get_sub_quads(self)
%GET_SUB_QUADS split the quad into subquads according to the reference points
%   NOTE: due to the needed ordering of the nodes for the matlab patch
%   function the node sequece of this function is different than usual. the
%   subquads are only used for the plotfunction of the fem model object.

% preallocate the array
sub_quads = zeros(self.order^2, 4);

% get the reference nodes ids
[~,sub_ids] = self.get_ref_points(self.order);

% loop over all the nodes of each element twice
sub_counter = 1;
for i=1:self.order
    for j = 1:self.order
        sub_quads(sub_counter,:) = [sub_ids(i,j),     sub_ids(i+1,j),...
                                    sub_ids(i+1,j+1), sub_ids(i,j+1)];
        sub_counter = sub_counter + 1;
    end
end
end

