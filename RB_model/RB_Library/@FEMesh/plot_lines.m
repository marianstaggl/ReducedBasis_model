function plot_lines(self)
%PLOT_LINES Plot the lines of the current setup
% get the unique ids of the lines
ids = unique([self.lines.id]);

% loop over all the ids and plot the lines
hold on;
for i=1:length(ids)
    
    % get the selector for the lines
    sel = [self.lines.id] == ids(i);
    
    % select the lines
    lines = self.lines(sel);

    % stich the lines together
    all_nodes = [{lines(:).pos}', repmat({[nan, nan]}, numel(lines), 1)];
    full_lines = cell2mat(reshape(all_nodes', [], 1));
    plot(full_lines(:,1), full_lines(:,2),'LineWidth', 2);
end

% plot the legend
legend(strsplit(num2str(ids)))
end
