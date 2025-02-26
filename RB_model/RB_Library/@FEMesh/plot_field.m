function plot_field( self, state_vec, vars, varargin )
%PLOT_FIELD plot a given state vector according to the own mesh
%   this function can be used to visualize the results of a simulation or
%   also certain lines of the discretized operator. which is also helpful a
%   lot of times

% parse the inputs
ip = inputParser();
ip.addParameter('marked_pts', false(size(state_vec)), @islogical);
ip.addParameter('tiled_layout', []);
ip.addParameter('link_colorbar', false);
ip.addParameter('parent', []);
ip.parse(varargin{:});

% clear the figure
if isempty(ip.Results.parent), c_fig = gcf; clf(c_fig);
else, c_fig = ip.Results.parent; end

% get the reference element
num_var = numel(vars);

% check if the nodes are in dl_arrays and extract if necessary
if isa(self.nodes, 'dlarray'), nodes = extractdata(self.nodes);
else, nodes = self.nodes; end
if isa(state_vec, 'dlarray'), state_vec = extractdata(state_vec); end

% reshape the state vector
state_vec = reshape(state_vec,[],1); 
state_vec = reshape(state_vec,[],num_var);
pts_vec = reshape(ip.Results.marked_pts, [], 1); 
pts_vec = reshape(pts_vec, [], num_var);

% preallocate the arrays
ref_subs = self.ref_elem.get_sub_quads();
num_subs = size(ref_subs,1);
subs = zeros(self.n_elem * (self.ref_elem.order-1)^2, 4);

% loop over all of the elements
s_c = 0;
for i=1:self.n_elem
    % insert the subquads
    subs(s_c + 1:s_c + num_subs,:) = reshape(self.elems(i,ref_subs),[],4);
    s_c = s_c + num_subs;
end

% plot a field for every variable
if isempty(ip.Results.tiled_layout)
    n_col = floor(sqrt(num_var));
    n_row = ceil(num_var/n_col);
else
    n_row = ip.Results.tiled_layout(1);
    n_col = ip.Results.tiled_layout(2);
end
t = tiledlayout(c_fig, n_row, n_col, 'TileSpacing','loose',...
    'Padding', 'loose');

% create the axis limits for the plot
x_range = [min(nodes(:,1)), max(nodes(:,1))];
x_limits = [x_range(1) - 0.05*diff(x_range),...
    x_range(2) + 0.05*diff(x_range)];
y_range = [min(nodes(:,2)), max(nodes(:,2))];
y_limits = [y_range(1) - 0.05*diff(y_range),...
    y_range(2) + 0.05*diff(y_range)];

for i=1:num_var
    
    % create a subplot
    ax(i) = nexttile(t); title(ax(i), vars(i)); hold(ax(i), 'on');
    patch(ax(i), 'Faces',[subs, subs(:,1)], 'Vertices', nodes,...
        'FaceVertexCData', state_vec(:,i),...
        'FaceColor','interp', 'LineStyle','none');

    % plot marked points if given
    if any(pts_vec(:,i))
        plot(ax(i), nodes(pts_vec(:,i), 1),...
            nodes(pts_vec(:,i), 2), 'ro'); 
    end
    
    % format the plot
    colorbar(ax(i)); 
    axis(ax(i), 'off'); hold(ax(i), 'off'); % axis(ax(i), 'equal', 'off'); 
    ax(i).set('XLim', x_limits, 'YLim', y_limits);

    if ip.Results.link_colorbar
        clim(ax(i), [min(state_vec, [], 'all'),...
            max(state_vec, [], 'all')]);
    end
end

% link the axes of the subplot
linkaxes(ax)

end

