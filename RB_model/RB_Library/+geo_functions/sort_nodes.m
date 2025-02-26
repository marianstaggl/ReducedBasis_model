function geo_lines = sort_nodes(geo_lines, resamp, patch_te)
%SORT_NODES sort the nodes of the lines

% remove all double points from the lines
for i=1:numel(geo_lines)
    geo_lines(i).pos = unique(geo_lines(i).pos, 'rows'); 
end

% loop over the blade lines
l_n = {geo_lines.name}; sel = find(contains(l_n, 'blade'));
for i=1:length(sel)
    
    % get a shortcut to the lines nodes
    pos = geo_lines(sel(i)).pos;
    
    % remove double points and use only southside of the blade
    pos = unique(pos, 'rows'); pos = pos(pos(:,2) <= 0, :);
    
    % sort the southside along the x axis
    [~,n_id] = sort(pos(:,1)); pos = pos(n_id, :);
    
    % shift last points in symmetry axis
    pos(1,2) = 0; pos(end,2) = 0;
    
    % do a resampling if needed
    if resamp, pos = resamp_line(pos, 200); end
    
    % put them together to close the loop
    geo_lines(sel(i)).pos = flipud([pos; flipud([pos(1:end-1,1),...
        -pos(1:end-1,2), pos(1:end-1,3)])]);
end

% patch the trailing edges if needed
if patch_te && resamp, geo_lines = patch_trail(geo_lines, 200, 0.11); end
end

function geo_lines = patch_trail(geo_lines, n_tr, pct)

% get a selector for the blades and tr nodes
b_id = find(contains({geo_lines.name}, 'blade'));
d_tr = round(n_tr*pct); tr_id = n_tr-d_tr:1:n_tr+d_tr;

% identify the middle contour
[tmp, tid] = sort([geo_lines(b_id).dist]); 
m_id = round(length(tmp)/2); b_m = geo_lines(b_id(tid(m_id)));

% get and normalize the median trailing edge
trm_x =  b_m.pos(tr_id,3); trm_y =  b_m.pos(tr_id,2);
trm_x = (trm_x - trm_x(1))/(max(trm_x) - min(trm_x)); 
trm_y = (trm_y - trm_y(1))/(abs(trm_y(end) - trm_y(1)));

% loop over the blades and replace the trailing edge
for i=1:length(b_id)
    % get the x and y values of the current tr_edge
    tr_x = geo_lines(b_id(i)).pos(tr_id,3);
    tr_y = geo_lines(b_id(i)).pos(tr_id,2);
    
    % get the x and y range for the tr edge
    dx = max(tr_x) - min(tr_x); dy = abs(tr_y(end) - tr_y(1));
    
    % dimensionalize tr_mean
    geo_lines(b_id(i)).pos(tr_id,3) = trm_x*dx + tr_x(1);
    geo_lines(b_id(i)).pos(tr_id,2) = trm_y*dy + tr_y(1);
end
end

function pos_new = resamp_line(pos, n)
% get a s vector for the contour line
s_old = cumsum(vecnorm(pos(2:end,:) - pos(1:end-1,:),2,2));
s_old = [0; s_old/s_old(end)]; s_new = functions.tanh_dist(2,n);

% interpolate x y and z values
pos_new(:,1) = reshape(makima(s_old, pos(:,1), s_new),[],1);
pos_new(:,2) = reshape(makima(s_old, pos(:,2), s_new),[],1);
pos_new(:,3) = reshape(makima(s_old, pos(:,3), s_new),[],1);
end

