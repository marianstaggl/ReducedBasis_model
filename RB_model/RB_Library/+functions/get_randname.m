function new_file = get_randname(folder)
%GET_RANDNAME generate a new name for a file

% loop until 1e10 to create new names
for i=1:1e10
    % generate a random hex number
    x = dec2hex(randi([1,1e10], 1));
    
    % set the datasetname
    new_file = fullfile(folder, ['dataset_', num2str(x), '.mat']);
    
    % check if the new name already exists
    if ~isfile(new_file), return; end
end

% if no new name can be found rise an error (this shouldn't happen)
error("couldn't find a new name...")
end

