function geo_lines = read_geofile(geo_file)
%READ_GEOFILES read the geometry files

% get the sheetnames of the file
sheets = sheetnames(geo_file);

% initialize the geo lines
geo_lines = struct('name', '', 'dist', nan, 'nodes', [], 'pos', []);
geo_lines = repmat(geo_lines, 1, numel(sheets));

% read the hub and shroud line
geo_lines(1).name = 'hub'; geo_lines(2).name = 'shroud';
geo_lines(1).pos = table2array(readtable(geo_file, 'Sheet', 'hub line'));
geo_lines(2).pos = table2array(readtable(geo_file, 'Sheet', 'shr line'));

% loop over the contour lines
sheets = sheets(contains(sheets, 'blade'));
for i=1:numel(sheets)
    % get the position of the line
    dist = strsplit(sheets{i},'0'); dist = str2double(dist{2});
    
    % save everything
    geo_lines(i+2).name = sheets{i}; geo_lines(i+2).dist = dist;
    geo_lines(i+2).pos = table2array(readtable(geo_file, 'Sheet', sheets{i}));
end
end

