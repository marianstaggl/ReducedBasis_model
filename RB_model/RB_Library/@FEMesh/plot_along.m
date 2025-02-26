function plot_along(self, line_id, state_vec, vars)
%PLOT_ALONG plot a profile along a specific line

% clear the figure
clf;

% get the number of variables
num_var = numel(vars);

% reshape the state vector
state_vec = reshape(state_vec,[],1); 
state_vec = reshape(state_vec,[],num_var);

% sort the nodes along the line
plt_vals = state_vec(self.lines(line_id).nodes,:);
plt_pos = self.get_s_vec(self.lines(line_id));
plt_pos = plt_pos/plt_pos(end);

% plot a field for every variable
for i=1:num_var
    % create a subplot
    ax(i) = subplot(num_var,1,i); title(vars(i))
    plot(plt_vals(:,i), plt_pos); grid on;
end

% link the axes of the subplot
% linkaxes(ax)

end

