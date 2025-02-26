function set_split_lines(self, threshold)
%SET_SPLIT_LINES split the lines by angle

% get the split lines
self.lines = self.static_split_lines(self.lines, threshold);
end

