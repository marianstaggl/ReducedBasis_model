function updateParameterValue(self)

selectedParam = self.ParamDropdown.Value;
self.Parameters.(selectedParam) = self.ParamInput.Value; % Update parameter value
self.updateTableData(); % Refresh table to reflect changes
end