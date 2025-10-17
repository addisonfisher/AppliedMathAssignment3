function generate_errors()
    
    t_ref = 0.492; 
    t_end = 10;

    [p_FE1_local, p_EM1_local, p_analytical1] = local_truncation(@rate_func01, t_ref, 1);
    [p_FE2_local, p_EM2_local, p_analytical2] = local_truncation(@rate_func02, t_ref, 2);
 
    [p_FE1_global_h, p_EM1_global_h, p_FE1_global_eval, p_EM1_global_eval] = ...
        global_truncation_error(@rate_func01, t_end, 1);

    [p_FE2_global_h, p_EM2_global_h, p_FE2_global_eval, p_EM2_global_eval] = ...
        global_truncation_error(@rate_func02, t_end, 2);
    
    fprintf('\n\n--- Local Truncation Error Scaling (p-values) ---\n');
    Methods_local = {'|X(t+h)-X(t)|'; 'Forward Euler'; 'Explicit Midpoint'};
    p_Test1_local = [p_analytical1; p_FE1_local; p_EM1_local];
    p_Test2_local = [p_analytical2; p_FE2_local; p_EM2_local];
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

