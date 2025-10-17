% To use, simply save this entire block of code as a single .m file
% and press the "Run" button in MATLAB.

function run_global_error_analysis()
    % set up the experiment parameters
    t_end1 = 10;   % end time for test 1
    t_end2 = 1000; % longer end time for test 2 to avoid machine precision issues

    % --- calculate p-values for test function 1 ---
    [p_fe1, p_be1, p_em1, p_im1] = calculate_p_values_for_test(...
        @rate_func01, @solution01, t_end1, 1, ...
        @backward_euler_step_analytical_rate01, @implicit_midpoint_step_analytical_rate01);

    % --- calculate p-values for test function 2 ---
    [p_fe2, p_be2, p_em2, p_im2] = calculate_p_values_for_test(...
        @rate_func02, @solution02, t_end2, [1; 0], ...
        @backward_euler_step_analytical_rate02, @implicit_midpoint_step_analytical_rate02);

    % --- display the results in a table ---
    fprintf('\n--- Global Error vs. Step Size (p-values) ---\n');
    methods = {'Forward Euler'; 'Backward Euler'; 'Explicit Midpoint'; 'Implicit Midpoint'};
    p_test1 = [p_fe1; p_be1; p_em1; p_im1];
    p_test2 = [p_fe2; p_be2; p_em2; p_im2];
    results_table = table(methods, p_test1, p_test2, 'VariableNames', {'Method', 'Estimate for p (test 1)', 'Estimate for p (test 2)'});
    disp(results_table);
end

% --- Main Calculation Function ---
function [p_fe, p_be, p_em, p_im] = calculate_p_values_for_test(rate_func, solution_func, t_end, x0, fast_be_step, fast_im_step)
    h_values = logspace(-2, 0, 25);
    errors = zeros(4, length(h_values)); % 4 methods
    x_analytical = solution_func(t_end);
    
    step_funcs = {@forward_euler_step, fast_be_step, @explicit_midpoint_step, fast_im_step};

    for i = 1:length(h_values)
        h = h_values(i);
        for j = 1:length(step_funcs)
            [~, x_numerical, ~, ~] = fixed_step_integration(rate_func, step_funcs{j}, [0, t_end], x0, h);
            errors(j, i) = norm(x_numerical(end, :)' - x_analytical);
        end
    end
    
    [p_fe, ~] = loglog_fit(h_values, errors(1, :));
    [p_be, ~] = loglog_fit(h_values, errors(2, :));
    [p_em, ~] = loglog_fit(h_values, errors(3, :));
    [p_im, ~] = loglog_fit(h_values, errors(4, :));
end

% --- Core Integration and Step Functions ---

function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    dx = rate_func_in(t, XA);
    XB = XA + dx * h;
    num_evals = 1;
end

function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    num_evals = 2;
    slope_initial = rate_func_in(t,XA);
    X_mid = XA + (h / 2) * slope_initial;
    t_mid = t + h / 2;
    slope_midpoint = rate_func_in(t_mid, X_mid);
    XB = XA + h * slope_midpoint;
end

function [t_list,X_list,h_avg, num_evals] = fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    t_start = tspan(1);
    t_end = tspan(2);
    num_steps = ceil((t_end - t_start) / h_ref);
    h_avg = (t_end - t_start) / num_steps;
    t_list = zeros(num_steps + 1, 1);
    X_list = zeros(num_steps + 1, size(X0, 1)); 
    t_list(1) = t_start;
    X_list(1, :) = X0'; 
    total_num_evals = 0;

    for i = 1:num_steps
        current_t = t_list(i);
        current_X = X_list(i, :)'; 
        [next_X, step_evals] = step_func(rate_func_in, current_t, current_X, h_avg);
        total_num_evals = total_num_evals + step_evals;
        t_list(i + 1) = t_list(i) + h_avg;
        X_list(i + 1, :) = next_X';
    end
    num_evals = total_num_evals;
end

% --- Rate and Solution Functions ---

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end

function X = solution02(t)
    X = [cos(t);sin(t)];
end

% --- Analytical Helpers for Implicit Methods (for speed) ---

% helpers for rate_func01
function [g_val, g_jac] = g_be_rate01(X_next, XA, h, t)
    rate_val = -5*X_next + 5*cos(t+h) - sin(t+h);
    g_val = XA + h * rate_val - X_next;
    g_jac = h * (-5) - 1;
end

function [g_val, g_jac] = g_im_rate01(X_next, XA, h, t)
    midpoint_X = 0.5 * (XA + X_next);
    midpoint_t = t + h/2;
    rate_val = -5*midpoint_X + 5*cos(midpoint_t) - sin(midpoint_t);
    g_val = XA + h * rate_val - X_next;
    g_jac = h * (-5) * 0.5 - 1;
end

function [XB,num_evals] = backward_euler_step_analytical_rate01(t, XA, h)
    g_with_jac = @(X_next_guess) g_be_rate01(X_next_guess, XA, h, t);
    solver_params = struct('dxmin',1e-10,'ftol',1e-10,'dxmax',1e8,'maxiter',50,'approx',0);
    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end

function [XB,num_evals] = implicit_midpoint_step_analytical_rate01(t, XA, h)
    g_with_jac = @(X_next_guess) g_im_rate01(X_next_guess, XA, h, t);
    solver_params = struct('dxmin',1e-10,'ftol',1e-10,'dxmax',1e8,'maxiter',50,'approx',0);
    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end

% helpers for rate_func02
function [g_val, g_jac] = g_be_rate02(X_next, XA, h)
    A = [0, -1; 1, 0];
    I = eye(2);
    g_val = XA + h * A * X_next - X_next;
    g_jac = h * A - I;
end

function [g_val, g_jac] = g_im_rate02(X_next, XA, h)
    A = [0, -1; 1, 0];
    I = eye(2);
    g_val = XA + h * A * (0.5 * (XA + X_next)) - X_next;
    g_jac = h * A * 0.5 - I;
end

function [XB,num_evals] = backward_euler_step_analytical_rate02(t, XA, h)
    g_with_jac = @(X_next_guess) g_be_rate02(X_next_guess, XA, h);
    solver_params = struct('dxmin',1e-10,'ftol',1e-10,'dxmax',1e8,'maxiter',50,'approx',0);
    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end

function [XB,num_evals] = implicit_midpoint_step_analytical_rate02(t, XA, h)
    g_with_jac = @(X_next_guess) g_im_rate02(X_next_guess, XA, h);
    solver_params = struct('dxmin',1e-10,'ftol',1e-10,'dxmax',1e8,'maxiter',50,'approx',0);
    [XB, ~, num_evals] = multi_newton_solver(g_with_jac, XA, solver_params);
end

% --- Utility Functions ---

function [p,k] = loglog_fit(x_regression,y_regression,varargin)
    if size(x_regression,1)==1, x_regression = abs(x_regression)'; end
    if size(y_regression,1)==1, y_regression = abs(y_regression)'; end
    Y = log(y_regression);
    X1 = log(x_regression);
    X2 = ones(length(X1),1);
    coeff_vec = regress(Y,[X1,X2]);
    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end

function [x, ~, num_evals] = multi_newton_solver(fun,x_guess,solver_params)
    num_evals = 0;
    dxmin = solver_params.dxmin;
    ftol = solver_params.ftol;
    dxmax = solver_params.dxmax;
    max_iter = solver_params.maxiter;
    approx = solver_params.approx;
    X = x_guess;
    delta_x = 1; 
    count = 0;

    if approx
        fval = fun(X);
        [J, jacobian_evals] = approximate_jacobian(fun, X);
        num_evals = num_evals + 1 + jacobian_evals;
    else
        [fval,J] = fun(X);
        num_evals = num_evals + 1;
    end

    while count < max_iter && norm(delta_x) > dxmin && norm(fval) > ftol && norm(delta_x) < dxmax
        count = count + 1;
        delta_x = -J\fval;
        X = X + delta_x;
        if approx
            fval = fun(X);
            [J, jacobian_evals] = approximate_jacobian(fun, X);
            num_evals = num_evals + 1 + jacobian_evals;
        else
            [fval,J] = fun(X);
            num_evals = num_evals + 1;
        end
    end
    x = X;
end

function [J, num_evals] = approximate_jacobian(fun,x)
    f_val = fun(x);
    num_inputs = length(x);
    num_outputs = length(f_val);
    J = zeros(num_outputs, num_inputs);
    h = 1e-6;
    for i = 1:num_inputs
        e = zeros(num_inputs, 1);
        e(i) = 1;
        J(:, i) = (fun(x + h*e) - fun(x - h*e)) / (2*h);
    end
    num_evals = 2 * num_inputs;
end