function val_data = validate_results(result_set,validation_set)
%VALIDATE_RESULTS Validate the results of a model with a dataset

% get the names fo the result and validation sets
if ~isequal({result_set.name}, {validation_set.name})
    error('result and validation set do not match...')
end

% loop over the validation set the total pressure
for i=1:numel(result_set)
    result_set(i).p_tot = get_ptot(result_set(i));
    validation_set(i).p_tot = get_ptot(validation_set(i));
end

% compare the velocities, pressure and density
vars = {'u', 'v', 'p', 'rho', 'p_tot'};

% loop over all results
for i=1:numel(vars)
    % collect the deviation in an array
    temp_val = [validation_set.(vars{i})];
    temp_test = [test_set.(vars{i})];
    
    % get the standard-deviation
    val_data.([vars{i}, '_std']) = mean(abs(temp_test - temp_val),'all');
    
    % get the max and min deviation
    val_data.([vars{i}, '_maxd']) = max(temp_test - temp_val, 'all');
    val_data.([vars{i}, '_mind']) = min(temp_test - temp_val, 'all');
    
    % get the overall correlation
    val_data.([vars{i}, '_corr']) = corr(...
        reshape(temp_val,[],1), reshape(temp_test,[],1));
end
end

function p_tot = get_ptot(res)
% calculate the total pressure at the nodes
p_tot = res.p + (res.rho/2)*(res.u.^2 + res.v.^2);
end
