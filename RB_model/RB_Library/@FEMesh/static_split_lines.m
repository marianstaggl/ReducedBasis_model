function n_lin = static_split_lines(lines, threshold)
%SPLIT_LINES split a line if the angle of two edges exceeds a value

% preallocate a new lines cell
n_lin = cell(1,numel(lines));
id = 0;

% loop over all of the lines
for i=1:numel(lines)
    
    % get the normalized vectors pescribing the edges (and their norms)
    vec = normc((lines(i).pos(1:end-1,:) - lines(i).pos(2:end,:))')';
    
    % calculate the angle between adjacent edges
    ang_m = acosd(sum(vec(1:end-1,:).*vec(2:end,:),2));
    
    % get the points where to split the edges
    sel = [1; find(ang_m > threshold)+1; length(lines(i).nodes)];
    
    % preallocate the strucarray
    n_lin{i} = repmat(struct('id', 0, 'nodes', [], 'pos', []),1,length(sel)-1);
    
    % loop over the split line parts and save them
    for j=1:length(sel)-1
        n_lin{i}(j).id = id; id = id + 1;
        n_lin{i}(j).pos = lines(i).pos(sel(j):sel(j+1),:);
        n_lin{i}(j).nodes = lines(i).nodes(sel(j):sel(j+1));
        % n_lin{i}(j).edges = lines(i).edges(sel(j):sel(j+1)-1);
        % n_lin{i}(j).elems = lines(i).elems(sel(j):sel(j+1)-1);
        % n_lin{i}(j).sub_id = lines(i).sub_id(:,sel(j):sel(j+1)-1);
        % n_lin{i}(j).n1 = lines(i).n1(sel(j):sel(j+1)-1,:);
        % n_lin{i}(j).n2 = lines(i).n2(sel(j):sel(j+1)-1,:);
    end
end

% reshape the cellarray to a strucarray
n_lin = cell2mat(n_lin);
end

