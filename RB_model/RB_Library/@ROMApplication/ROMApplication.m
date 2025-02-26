classdef ROMApplication < handle
    properties
        % Define UI components
        Figure
        PlotPanel
        Axes
        RO_Model
        ROM_Function
        ParamDropdown
        ParamInput
        PlotButton
        ParamTable

        % Properties for parameter values
        Parameters = struct();
    end
    
    methods
        % Constructor to initialize the GUI components
        function self = ROMApplication(ro_model, run_fun, par_list)
            % store the ro_model and the parameter list
            self.RO_Model = ro_model;
            self.ROM_Function = run_fun;
            par_val = num2cell(zeros(1, numel(par_list)));
            temp = reshape(cat(1, par_list, par_val), [], 1);
            self.Parameters = struct(temp{:});
            
            % Create the main figure
            self.Figure = uifigure('Name', 'ROMApplication');
            
            % Create a panel for plotting, taking most of the figure width
            self.PlotPanel = uipanel(self.Figure, 'Position', [20, 20, 400, 380]);
            c_ax = axes(self.PlotPanel); % create an axes to achieve desired resize behaviour
            axis(c_ax, 'off');
            
            % Create a dropdown menu for parameter selection on the right
            self.ParamDropdown = uidropdown(self.Figure, ...
                'Position', [440, 375, 100, 22], ...
                'Items', par_list, ...
                'Value', par_list{1}, ...
                'ValueChangedFcn', @(src, ~) self.updateInputField(src.Value));
            
            % Create an input field for parameter values below the dropdown
            self.ParamInput = uieditfield(self.Figure, 'numeric', ...
                'Position', [440, 340, 100, 22], ...
                'Value', self.Parameters.(par_list{1}), ...
                'ValueChangedFcn', @(~, ~) self.updateParameterValue());
            
            % Create a button to execute the function under the input field
            self.PlotButton = uibutton(self.Figure, 'push', ...
                'Position', [440, 305, 100, 22], ...
                'Text', 'Plot', ...
                'ButtonPushedFcn', @(~, ~) self.plotFunction());
            
            % Create a table to display parameter values below the plot button
            self.ParamTable = uitable(self.Figure, ...
                'Position', [440, 20, 100, 140], ...
                'ColumnName', {'Parameter', 'Value'}, ...
                'Data', self.getParameterData());
        end
    end
end
