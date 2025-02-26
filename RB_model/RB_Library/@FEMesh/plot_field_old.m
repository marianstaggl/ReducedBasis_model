function plot_field( self, state_vec, vars, varargin )
%PLOT_FIELD plot a given state vector according to the own mesh
%   this function can be used to visualize the results of a simulation or
%   also certain lines of the discretized operator. which is also helpful a
%   lot of times

% parse the inputs
ip = inputParser();
ip.addParameter('marked_pts', false(size(state_vec)), @islogical);
ip.parse(varargin{:});

% clear the figure
clf;

% get the reference element
num_var = numel(vars);

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
n_col = floor(sqrt(num_var));
n_row = ceil(num_var/n_col); 
for i=1:num_var
    
    % create a subplot
    ax(i) = subplot(n_row,n_col,i); title(vars(i)); hold on;
    patch('Faces',[subs, subs(:,1)], 'Vertices', self.nodes,...
        'FaceVertexCData', state_vec(:,i),...
        'FaceColor','interp', 'LineStyle','none');

    % plot marked points if given
    if any(pts_vec(:,i))
        plot(self.nodes(pts_vec(:,i), 1),...
            self.nodes(pts_vec(:,i), 2), 'ro'); 
    end
    
    % format the plot
    colorbar(); axis equal; axis off; hold off;
end

% link the axes of the subplot
linkaxes(ax)

end

