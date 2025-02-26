function data = getParameterData(self)
params = fieldnames(self.Parameters);
values = struct2cell(self.Parameters);
data = [params, values];
end