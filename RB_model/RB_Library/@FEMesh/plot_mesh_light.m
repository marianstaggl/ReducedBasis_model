function plot_mesh_light(self, varargin)
%PLOT_MESH Plot the mesh of the current setup

% parse the inputs
ip = inputParser();
ip.addParameter('marked_elem', []);
ip.parse(varargin{:});

% get the corner nodes of the current element
[~, ids] = FElement.get_ref_points( self.ref_elem.order );
corners = [ids(1,1), ids(end,1), ids(end,end), ids(1,end), ids(1,1)];

% hold the figure
hold on; xlabel('x'); ylabel('y'); title('mesh topology'); axis equal;

% loop over all of the elements and plot them
for i=1:size(self.elems,1)
    % get the current elements points
    curr_el = self.elems(i,:);
    curr_n = self.nodes(curr_el,:);
    
    % plot the elements edges
    plot(curr_n(corners,1), curr_n(corners,2), 'Color', 'k');
end

% use patch to fill the marked elements
patch('Faces', self.elems(ip.Results.marked_elem, corners),...
    'Vertices', self.nodes,'FaceColor','red', 'FaceAlpha', 0.5,...
    'LineStyle','none');
end

