function calculate_global_error_scaling()
    % This script calculates the scaling of the global truncation error with step size
    % for several numerical methods and test functions. It then displays the results
    % in a table.

    % Add the necessary subdirectories to the path
    addpath('./Applied-Math-Assignment-2');

    % Define the test parameters
    t_end = 10;

    % --- Test Function 1 ---
    [p_FE1, p_BE1, p_EM1, p_IM1] = calculate_p_values(@rate_func01, @solution01, [0, t_end], 1);

    % --- Test Function 2 ---
    [p_FE2, p_BE2, p_EM2, p_IM2] = calculate_p_values(@rate_func02, @solution02, [0, t_end], [1; 0]);


    % --- Display Results in a Table ---
    fprintf('\n--- Global Error vs. Step Size (p-values) ---\n');
    Methods = {'Forward Euler'; 'Backward Euler'; 'Explicit Midpoint'; 'Implicit Midpoint'};
    p_Test1 = [p_FE1; p_BE1; p_EM1; p_IM1];
    p_Test2 = [p_FE2; p_BE2; p_EM2; p_IM2];
    results_table = table(Methods, p_Test1, p_Test2, 'VariableNames', {'Method', 'p (test 1)', 'p (test 2)'});
    disp(results_table);

end

function [p_FE, p_BE, p_EM, p_IM] = calculate_p_values(rate_func_in, solution_func, t_span, X0)
    % This helper function calculates the p-values for a given test case.

    % Range of reference step sizes to test
    h_ref_list = logspace(-3, 0, 50);

    % Allocate arrays to store results
    h_list_FE = zeros(size(h_ref_list));
    error_FE = zeros(size(h_ref_list));

    h_list_BE = zeros(size(h_ref_list));
    error_BE = zeros(size(h_ref_list));

    h_list_EM = zeros(size(h_ref_list));
    error_EM = zeros(size(h_ref_list));

    h_list_IM = zeros(size(h_ref_list));
    error_IM = zeros(size(h_ref_list));


    X_analytical_final = solution_func(t_span(2));

    % Loop through the different step sizes
    for i = 1:length(h_ref_list)
        h_ref = h_ref_list(i);

        % --- Forward Euler ---
        [~, X_list_FE, h_avg_FE, ~] = fixed_step_integration(rate_func_in, @forward_euler_step, t_span, X0, h_ref);
        X_final_FE = X_list_FE(end, :)';
        error_FE(i) = norm(X_final_FE - X_analytical_final);
        h_list_FE(i) = h_avg_FE;

        % --- Backward Euler ---
        [~, X_list_BE, h_avg_BE, ~] = fixed_step_integration(rate_func_in, @backward_euler_step, t_span, X0, h_ref);
        X_final_BE = X_list_BE(end, :)';
        error_BE(i) = norm(X_final_BE - X_analytical_final);
        h_list_BE(i) = h_avg_BE;

        % --- Explicit Midpoint ---
        [~, X_list_EM, h_avg_EM, ~] = fixed_step_integration(rate_func_in, @explicit_midpoint_step, t_span, X0, h_ref);
        X_final_EM = X_list_EM(end, :)';
        error_EM(i) = norm(X_final_EM - X_analytical_final);
        h_list_EM(i) = h_avg_EM;

        % --- Implicit Midpoint ---
        [~, X_list_IM, h_avg_IM, ~] = fixed_step_integration(rate_func_in, @implicit_midpoint_step, t_span, X0, h_ref);
        X_final_IM = X_list_IM(end, :)';
        error_IM(i) = norm(X_final_IM - X_analytical_final);
        h_list_IM(i) = h_avg_IM;
    end

    % Fit lines to the log-log plot to find the p-values
    [p_FE, ~] = loglog_fit(h_list_FE, error_FE);
    [p_BE, ~] = loglog_fit(h_list_BE, error_BE);
    [p_EM, ~] = loglog_fit(h_list_EM, error_EM);
    [p_IM, ~] = loglog_fit(h_list_IM, error_IM);

end