function new_mesh = copy(self)
%COPY Return a deep copy of the object

% get the current working directory
temp_file = [cd '\temp.mat'];

% save the cfd results and their parameter
save(temp_file, 'self','-v7.3','-nocompression')

% reload the written object
new_mesh = load(temp_file); new_mesh = new_mesh.self;

% clean up afterwards
delete(temp_file)

end

