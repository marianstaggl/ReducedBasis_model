function save_par(par, varargin)
%SAVE_PAR save the parametrization as excel file

% set the default names
c_names = repmat({'none'}, 1, par(1).lines);

% parse the inputs
ip = inputParser();
ip.addParameter('names', c_names)
ip.addParameter('filename', 'bc_data.xlsx')
ip.parse(varargin{:});

% put the parameters into tables
for i=1:numel(par)
   % convert the data into a table
   T = array2table(par(i).data,'VariableNames',par(i).columns);
   
   % add the names to the table
   T.names = reshape(ip.Results.names,[],1);
   
   % move the name to the front
   T = movevars(T,'names','Before',par(i).columns{1}); 
   
   % save the table as a excel sheet
   writetable(T, ip.Results.filename, 'Sheet', par(i).name);
end
end

