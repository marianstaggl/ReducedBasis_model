%% clear the workspace and define parameters
% clear the workspace and add the library
clc; clear all; close all; 
addpath('..\..\RB_Library')
root_path = '..\..\RB_data\rom_root_sudden_expansion_channel';

%% setup the reduced model
% initialize the reduced order model
ro_model = ROModel_isa(root_path);
ro_model.set_base_fun('em19_gm4');
ro_model.set_sys_mat();
ro_model.set_deim_elem();

% load the geometrical base functions
geo_base = load('+util_fun\geo_base.mat').geo_base;
geo_fun = @(x) reshape(geo_base(:,1) + x*geo_base(:,2), [], 2);

%% specify the target value (can be a scalar or a vector)
target_value = 8e-3 * ones(200, 1); 
target_value(1:50) = 0 * target_value(1:50);

%% evaluate the reduced model
% get the initial parameters
init_par = [0.8, 0.5, 0];

% define the objective function
eval_fun = @(ro_model, par_vec) util_fun.evaluate_model(ro_model, par_vec, geo_fun);
obj_fun = @(sol) mean((reshape(util_fun.obj_velprofile(ro_model, sol), [], 1)...
    - target_value).^2, "all");

% run the optimization loop
param = ro_model.run_optimization(eval_fun, obj_fun, init_par,...
    'n_iter', 5, 'extrap_penalty', 1);

%% plot the final result
% evaluate the parameter combination and plot the result
s = util_fun.evaluate_model(ro_model, param, geo_fun);
ro_model.fe_model.plot_field(extractdata(s));

% plot the result of the optimization
figure; hold on;
plot(extractdata(util_fun.obj_velprofile(ro_model, s)),linspace(0, 1, 200)); 
plot(target_value,linspace(0, 1, 200)); grid on;
legend({'optimization result', 'target value'})

opti_case.target_value = target_value;
opti_case.optim_results = extractdata(util_fun.obj_velprofile(ro_model, s));
save('opti_case', 'opti_case')