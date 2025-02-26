function plotFunction(self)
% run the reduced order model on the given parameters
new_sol = self.ROM_Function(...
    self.RO_Model, self.Parameters);

% Clear any previous plot
self.RO_Model.fe_model.plot_field(new_sol,...
    'parent', self.PlotPanel);
end