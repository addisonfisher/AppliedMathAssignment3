function generate_errors()
    
    t_ref = 0.492; 
    t_end = 10;

    [p_FE1_local, p_EM1_local, p_analytical1] = local_truncation(@rate_func01, t_ref, 1);
    [p_FE2_local, p_EM2_local, p_analytical2] = local_truncation(@rate_func02, t_ref, 2);
 
    [p_FE1_global_h, p_EM1_global_h, p_FE1_global_eval, p_EM1_global_eval] = ...
        global_truncation_error(@rate_func01, t_end, 1);

    [p_FE2_global_h, p_EM2_global_h, p_FE2_global_eval, p_EM2_global_eval] = ...
        global_truncation_error(@rate_func02, t_end, 2);
   

    % setup for local truncation error calculation
    h = logspace(-5, 1, 100);
    error_list_be1 = [];
    error_list_im1 = [];
    error_list_be2 = [];
    error_list_im2 = [];

    % calculate errors for backward euler and implicit midpoint
    for i = 1:length(h)
        % test function 1
        x_be1 = solution01(t_ref);
        x_h_be1 = solution01(t_ref + h(i));
        g_x_be1 = backward_euler_step(@rate_func01, t_ref, x_be1, h(i));
        error_list_be1(end + 1) = norm(g_x_be1 - x_h_be1);

        x_im1 = solution01(t_ref);
        x_h_im1 = solution01(t_ref + h(i));
        g_x_im1 = implicit_midpoint_step(@rate_func01, t_ref, x_im1, h(i));
        error_list_im1(end + 1) = norm(g_x_im1 - x_h_im1);

        % test function 2
        x_be2 = solution02(t_ref);
        x_h_be2 = solution02(t_ref + h(i));
        g_x_be2 = backward_euler_step(@rate_func02, t_ref, x_be2, h(i));
        error_list_be2(end + 1) = norm(g_x_be2 - x_h_be2);

        x_im2 = solution02(t_ref);
        x_h_im2 = solution02(t_ref + h(i));
        g_x_im2 = implicit_midpoint_step(@rate_func02, t_ref, x_im2, h(i));
        error_list_im2(end + 1) = norm(g_x_im2 - x_h_im2);
    end
    
    % fit data to get p-values
    [p_BE1_local, ~] = loglog_fit(h, error_list_be1);
    [p_IM1_local, ~] = loglog_fit(h, error_list_im1);
    [p_BE2_local, ~] = loglog_fit(h, error_list_be2);
    [p_IM2_local, ~] = loglog_fit(h, error_list_im2);
    
    % --- end of added code ---

    fprintf('\n\n--- Local Truncation Error Scaling (p-values) ---\n');
    Methods_local = {'|X(t+h)-X(t)|'; 'Forward Euler'; 'Backward Euler'; 'Explicit Midpoint'; 'Implicit Midpoint'};
    p_Test1_local = [p_analytical1; p_FE1_local; p_BE1_local; p_EM1_local; p_IM1_local];
    p_Test2_local = [p_analytical2; p_FE2_local; p_BE2_local; p_EM2_local; p_IM2_local];
    results_table_local = table(Methods_local, p_Test1_local, p_Test2_local, 'VariableNames', {'Method', 'Estimate for p (test 1)', 'Estimate for p (test 2)'});
    disp(results_table_local);

    fprintf('\n--- Global Error vs. Step Size (p-values) ---\n');
    Methods_global_h = {'Forward Euler'; 'Explicit Midpoint'};
    p_Test1_global_h = [p_FE1_global_h; p_EM1_global_h];
    p_Test2_global_h = [p_FE2_global_h; p_EM2_global_h];
    table_global_h = table(Methods_global_h, p_Test1_global_h, p_Test2_global_h, 'VariableNames', {'Method', 'p (test 1)', 'p (test 2)'});
    disp(table_global_h);

    fprintf('\n--- Global Error vs. Function Evals (p-values) ---\n');
    Methods_global_eval = {'Forward Euler'; 'Explicit Midpoint'};
    p_Test1_global_eval = [p_FE1_global_eval; p_EM1_global_eval];
    p_Test2_global_eval = [p_FE2_global_eval; p_EM2_global_eval];
    table_global_eval = table(Methods_global_eval, p_Test1_global_eval, p_Test2_global_eval, 'VariableNames', {'Method', 'p (test 1)', 'p (test 2)'});
    disp(table_global_eval);
end