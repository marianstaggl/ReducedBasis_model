function plot_boundarys(self)
%PLOT_BOUNDARYS Plot the boundarys of the current setup
%   to check the boundary conditions of the current setup they are plotted
%   here.

% get the default plot colors of matlab
cmap = [0.6350, 0.0780, 0.1840;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560;...
    0.4660, 0.6740, 0.1880;...
    0.3010, 0.7450, 0.9330];

% get the number of system variables
sys_vars = self.sys_variables; num_vars = numel(sys_vars);

% reshape the dir boundary conditions to account for different variables
d_sel = reshape(self.dir_bounds(:,1), [], num_vars);
d_val = reshape(self.dir_bounds(:,2), [], num_vars);

% reshape the neu boundary conditions to account for different variables
n_sel = reshape(self.neu_bounds(:,1), [], num_vars);
n_val = reshape(self.neu_bounds(:,2), [], num_vars);

% reshape the neu boundary conditions to account for different variables
r_sel = reshape(self.rob_bounds(:,1), [], num_vars);
r_val = reshape(self.rob_bounds(:,2), [], num_vars);

% loop over the variables and plot the dirichlet boundary conditions
for i=1:num_vars
    
    % get the plot indices
    pidx1 = (i-1)*3 + 1; pidx2 = (i-1)*3 + 2; pidx3 = (i-1)*3 + 3;
    
    % plot the dirichlet boundary conditions
    h = subplot(num_vars, 3, pidx1); cla; title([sys_vars{i} ': dirichlet boundarys']);
    plot_bounds(h, cmap, self.mesh.nodes, [d_sel(:,i), d_val(:,i)]);
    
    % plot the neumann boundary conditions
    h = subplot(num_vars, 3, pidx2); cla; title([sys_vars{i} ': neumann boundarys'])
    plot_bounds(h, cmap, self.mesh.nodes, [n_sel(:,i), n_val(:,i)]);
    
    % plot the neumann boundary conditions
    h = subplot(num_vars, 3, pidx3); cla; title([sys_vars{i} ': robin boundarys'])
    plot_bounds(h, cmap, self.mesh.nodes, [r_sel(:,i), r_val(:,i)]);
end


end

function plot_bounds(h, cmap,  nodes, bounds)

% plot the nodes of the mesh
hold on; plot(nodes(:,1), nodes(:,2),'.', 'Parent', h);

% get all the unique values of the dirichlet bounds
b_val = unique(bounds(:,2)); legnd = cell(1,0); lc = 1;
for i=1:numel(b_val)
    
    % select the boundary nodes
    sel = (bounds(:,2) == b_val(i)) & bounds(:,1);
    
    % add the value to the legend
    if any(sel)
        legnd{lc} = ['bound value: ' num2str(b_val(i))];
        lc = lc + 1;
    end
    
    % plot the boundary nodes
    plot(nodes(sel,1), nodes(sel,2), 'o', 'Color',...
        cmap(rem(i,size(cmap,1))+1, :), 'Parent', h);
end

% plot the legend of the boundary values
legend(horzcat({'grid pts'}, legnd), 'Location', 'SouthWest');

end

