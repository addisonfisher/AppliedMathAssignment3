function global_error_func()
    % set up the experiment parameters
    t_end = 100;
    addpath('Applied-Math-Assignment-2');
    % calculate p-values for the first test function
    [p_be1, p_im1] = calculate_for_test_function(@rate_func01, @solution01, t_end, 1);

    % calculate p-values for the second test function
    [p_be2, p_im2] = calculate_for_test_function(@rate_func02, @solution02, t_end, [1; 0]);

    % display the results in a table
    fprintf('\n--- Global Error vs. Function Evals (p-values) ---\n');
    methods = {'Backward Euler'; 'Implicit Midpoint'};
    p_test1 = [p_be1; p_im1];
    p_test2 = [p_be2; p_im2];
    results_table = table(methods, p_test1, p_test2, 'VariableNames', {'Method', 'Estimate for p (test 1)', 'Estimate for p (test 2)'});
    disp(results_table);
end

function [p_be, p_im] = calculate_for_test_function(rate_func, solution_func, t_end, x0)
    % set up a range of step sizes
    h_values = logspace(-3, 0, 25);
    
    % pre-allocate arrays
    be_errors = zeros(1, length(h_values));
    im_errors = zeros(1, length(h_values));
    be_evals = zeros(1, length(h_values));
    im_evals = zeros(1, length(h_values));
    
    % get the analytical solution at the end time
    x_analytical = solution_func(t_end);

    % create handles for the fast analytical step functions
    fast_be_step = @(rf, t, xa, h) backward_euler_step_analytical(t, xa, h);
    fast_im_step = @(rf, t, xa, h) implicit_midpoint_step_analytical(t, xa, h);

    % loop through each step size
    for i = 1:length(h_values)
        h = h_values(i);
        
        % backward euler
        [~, x_be, ~, num_evals_be] = fixed_step_integration(rate_func, fast_be_step, [0, t_end], x0, h);
        be_errors(i) = norm(x_be(end, :)' - x_analytical);
        be_evals(i) = num_evals_be;
        
        % implicit midpoint
        [~, x_im, ~, num_evals_im] = fixed_step_integration(rate_func, fast_im_step, [0, t_end], x0, h);
        im_errors(i) = norm(x_im(end, :)' - x_analytical);
        im_evals(i) = num_evals_im;
    end
    
    % perform log-log fit to find p-values
    [p_be, ~] = loglog_fit(be_evals, be_errors);
    [p_im, ~] = loglog_fit(im_evals, im_errors);
end


%%this uses helper functions due to speed issues with other implementations
%had an original function written but used some chatgpt to help optimize so
%it actually converges.

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