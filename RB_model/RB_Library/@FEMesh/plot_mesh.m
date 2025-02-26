function plot_mesh(self)
%PLOT_MESH Plot the mesh of the current setup

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
    
    % plot the elements points
    plot(curr_n(:,1), curr_n(:,2), ':o')
    
    % plot the elements edges
    plot(curr_n(corners,1), curr_n(corners,2), 'Color', 'k');
end
end

