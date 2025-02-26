function ids = static_select_lines(lines, cid)
%SELECT_LINES use ginput to select a subset of the given lines

while true
    % clear the figure and set it to hold
    clf; hold on; grid on; axis equal;
    
    % loop over the lines and plot them
    plot_lines(lines, cid, '-');
    
    % use ginput to select the lines via userinput
    [x,y] = functions.ginput_zoom();
    
    % find the closest points for every line (use unique)
    ids = unique(closest_lines(lines, cid, [x,y]));
    
    % plot the selected lines with dots to check result
    plot_lines(lines(ids), cid, 'o'); pause();
    
    % check if we want to break the loop
    br = input('reject the selection? [y/n]', "s");
    
    % check if everything is fine
    if isequal('n',lower(br)), break; end
end
end

function ids = closest_lines(lines, cid, p)
% preallocate the arrays
min_d = zeros(1,numel(lines)); ids = zeros(1,size(p,1));

% loop over the lines and get the distance
for i=1:size(p,1)
    % get the smallest dist to the lines
    for j=1:numel(lines), min_d(j) = min(vecnorm(lines(j).pos(:,cid) - p(i,:),2,2)); end
    
    % save the id of the smallest min_d value
    [~, ids(i)] = min(min_d);
end
end

function plot_lines(lines, cid, mar)
% loop over the lines and plot them
for i=1:numel(lines)
    plot(lines(i).pos(:,cid(1)), lines(i).pos(:,cid(2)), mar);
end
end
