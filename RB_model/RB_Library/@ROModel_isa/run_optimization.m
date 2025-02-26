function parameter = run_optimization(self, eval_fun, obj_fun, init_param, varargin)

% parse the input parameters
ip = inputParser();
ip.addParameter('n_iter', 15);
ip.addParameter('conv_criterion', 1e-12);
ip.addParameter('extrap_penalty', 1);
ip.parse(varargin{:});

% specify a loss function for the extrapolation penalty
extrap_loss = @(param) ip.Results.extrap_penalty * ...
    (sum(relu(- param) + relu(param - 1)));

% create a wrapper for the loss function
target_fun = @(x) obj_fun(eval_fun(self, x)) + extrap_loss(x);
loss_fun = @(x) dlfeval(@model_loss, target_fun, x);

% run the optimizer for the given loss function
conv_criterion = ip.Results.conv_criterion;
parameter = dlarray(init_param);
solver_state = lbfgsState('HistorySize',50);
for i=1:ip.Results.n_iter
    % initialize the solver
    [parameter, solver_state] = lbfgsupdate(parameter,...
        loss_fun, solver_state);

    current_loss = extractdata(solver_state.Loss);
    disp(['iteration ', num2str(i), ' out of ', num2str(ip.Results.n_iter),...
        ' | current loss: ', num2str(current_loss)])

    % get the convergence of the loss
    loss_convergence(i) = current_loss;
    loss_diff = diff(loss_convergence);

    plot_convergence(loss_convergence)

    conv_check = current_loss < conv_criterion & ...
        min(abs(loss_diff)) < 0.1 * conv_criterion;
    if conv_check, break; end
end
end

function [loss, gradients] = model_loss(loss_fun, param)
loss = loss_fun(param);
gradients = dlgradient(loss, param);
end

function plot_convergence(loss_values)
clf; semilogy(loss_values);
grid on; title('loss convergence')
xlabel('# of iterations');
ylabel('loss value'); drawnow;
end