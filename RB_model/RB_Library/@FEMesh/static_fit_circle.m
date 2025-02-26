function [xm, ym, r] = static_fit_circle(lines, cid)
%FIT_CIRCLE use ginput to fit a circle into a line

% get all the nodes into one array
all_nodes = vertcat(lines.pos);

% loop until the result is ok
while true
    % clear the figure and set it to hold
    clf; hold on; grid on; axis equal;
    
    % loop over the lines and plot them
    plot_lines(lines, cid, '-');
    
    % use ginput to select the lines via userinput
    [x,y] = functions.ginput_zoom();
    
    % find the closest points for every line (use unique)
    [xcl, ycl] = closest_point(all_nodes(:,cid), [x,y]);
    
    % check if there are at least 3 points
    if numel(xcl) < 3, continue; end
    
    % fit a circle to the selected points and plot
    [xm, ym, r] = functions.circle_fit(xcl, ycl); viscircles([xm, ym], r)
    
    % check if we want to break the loop
    br = input('reject the selection? [y/n]', "s");
    
    % check if everything is fine
    if isequal('n',lower(br)), break; end
end
end

function [x,y] = closest_point(xy_nodes, sel_nodes)
% get the closest nodes to the selected ones

for i=1:size(sel_nodes,1)
    % get the distances to the nodes
    d_vec = vecnorm(xy_nodes-sel_nodes(i,:),2,2);
    
    % get the closest node
    [~,id] = min(d_vec); 
    x(i) = xy_nodes(id,1); y(i) = xy_nodes(id,2);
end
end

function plot_lines(lines, cid, mar)
% loop over the lines and plot them
for i=1:numel(lines)
    plot(lines(i).pos(:,cid(1)), lines(i).pos(:,cid(2)), mar);
end
end