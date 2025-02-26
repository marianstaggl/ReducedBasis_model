function [par, names] = load_par(filename)
%LOAD_PAR load the parameters for bcs

% get the sheet names
sheets = sheetnames(filename);

% preallocate the parobjects
par(numel(sheets)) = DataObj();

% loop over the sheet names
for i=1:numel(sheets)
    % read the excel sheet
    T = readtable(filename,'Sheet', sheets(i));
    
    % set the names return
    names = reshape(T.names,1,[]);
    data = table2array(T(:,~ismember(T.Properties.VariableNames, 'names')));
    coln = T.Properties.VariableNames(~ismember(T.Properties.VariableNames, 'names'));
    
    % set the name of the parameter
    par(i).name = char(sheets(i));
    par(i).add_column(coln, data)
end
end

