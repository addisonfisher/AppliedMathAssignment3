function plot_global_truncation_all_methods()
    % set up parameters for the experiment
    t_span = [0, 10];
    x0 = 1;
    % you can now use more points without a long wait
    h_values = logspace(-3, 0, 25); 

    % create handles to the fast, analytical versions of the implicit step functions
    fast_be_step = @(rate_func_in, t, XA, h) backward_euler_step_analytical(t, XA, h);
    fast_im_step = @(rate_func_in, t, XA, h) implicit_midpoint_step_analytical(t, XA, h);

    % define the methods to be plotted, using the fast versions for implicit methods
    step_funcs = {@forward_euler_step, fast_be_step, @explicit_midpoint_step, fast_im_step};
    method_names = {'Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint'};
    num_methods = length(step_funcs);
    
    % pre-allocate the errors matrix for speed
    errors = zeros(num_methods, length(h_values));

    % analytical solution at the end point
    x_analytical = solution01(t_span(2));

    % calculate errors for each method and h value
    for i = 1:num_methods
        for j = 1:length(h_values)
            h = h_values(j);
            % run the integration
            [~, x_numerical, ~, ~] = fixed_step_integration(@rate_func01, step_funcs{i}, t_span, x0, h);
            % calculate the error at the final time step
            errors(i, j) = abs(x_numerical(end) - x_analytical);
        end
    end

    % create the log-log plot
    figure;
    hold on;
    
    % plot data points for each method
    loglog(h_values, errors(1, :), '-', 'DisplayName', method_names{1});
    loglog(h_values, errors(2, :), '-', 'DisplayName', method_names{2});
    loglog(h_values, errors(3, :), '-', 'DisplayName', method_names{3});
    loglog(h_values, errors(4, :), '-', 'DisplayName', method_names{4});
    set(gca, 'XScale', 'log', 'YScale', 'log')

    % format the plot
    title('Global Truncation Error vs. Step Size');
    xlabel('Step Size (h)');
    ylabel('Global Truncation Error');
    legend('show', 'Location', 'southeast');
    hold off;
end


% --- local helper functions for speed ---

% helper for backward euler g-function and its analytical jacobian
function [g_val, g_jac] = g_be_rate01(X_next, XA, h, t)
    rate_val = -5*X_next + 5*cos(t+h) - sin(t+h);
    g_val = XA + h * rate_val - X_next;
    % analytical jacobian of g is h * j_f - i
    g_jac = h * (-5) - 1;
end

% helper for implicit midpoint g-function and its analytical jacobian
function [g_val, g_jac] = g_im_rate01(X_next, XA, h, t)
    midpoint_X = 0.5 * (XA + X_next);
    midpoint_t = t + h/2;
    rate_val = -5*midpoint_X + 5*cos(midpoint_t) - sin(midpoint_t);
    g_val = XA + h * rate_val - X_next;
    % analytical jacobian of g is h * j_f * 0.5 - i
    g_jac = h * (-5) * 0.5 - 1;
end

% fast backward euler step that uses the analytical jacobian
function [XB,num_evals] = backward_euler_step_analytical(t, XA, h)
    g_with_jac = @(X_next_guess) g_be_rate01(X_next_guess, XA, h, t);
    
    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.maxiter = 50; % reduced maxiter
    solver_params.approx = 0; % use analytical jacobian

    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end

% fast implicit midpoint step that uses the analytical jacobian
function [XB,num_evals] = implicit_midpoint_step_analytical(t, XA, h)
    g_with_jac = @(X_next_guess) g_im_rate01(X_next_guess, XA, h, t);

    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.maxiter = 50; % reduced maxiter
    solver_params.approx = 0; % use analytical jacobian

    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end