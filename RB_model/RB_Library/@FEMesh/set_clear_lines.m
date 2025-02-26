function set_clear_lines(self, line_id)
%CLEAR_LINES remove intersections of the lines
%    sometimes the lines intersect at the corners. this function can be
%    used to remove those redundancies

% identify all the lines with the searched id
sel = ([self.lines.id] == line_id);
idx = 1:1:length(self.lines); idx = idx(~sel);

% identify all nodes which are presend in the lines
all_n = [self.lines(sel).nodes];

% loop over the indices of the remaining lines
for i=1:length(idx)
    subsel = ismember(self.lines(idx(i)).nodes, all_n);
    self.lines(idx(i)).nodes(subsel) = [];
end
end

