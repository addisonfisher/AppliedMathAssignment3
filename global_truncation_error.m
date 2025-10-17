function [p_FE_h, p_EM_h, p_FE_eval, p_EM_eval] = global_truncation_error(rate_func_in, t_end, flag)
    
    if flag == 1
        solution_func = @solution01;
        X0 = 1;
        t_span = [0, t_end];
    elseif flag == 2
        solution_func = @solution02;
        X0 = [1; 0];
        t_span = [0, t_end];
    end

    %range of reference step sizes to test
    h_ref_list = logspace(-3, 0, 50);

    %allocate arrays to store results
    h_list = zeros(size(h_ref_list));
    error_FE = zeros(size(h_ref_list));
    error_EM = zeros(size(h_ref_list));
    num_evals_FE = zeros(size(h_ref_list));
    num_evals_EM = zeros(size(h_ref_list));

    X_analytical_final = solution_func(t_span(2));

    %loop through
    for i = 1:length(h_ref_list)
        h_ref = h_ref_list(i);

        %full integration with forward ruler
        [~, X_list_FE, h_avg, evals_FE] = fixed_step_integration(rate_func_in, @forward_euler_step, t_span, X0, h_ref);
        X_final_FE = X_list_FE(end, :)';
        error_FE(i) = norm(X_final_FE - X_analytical_final);
        num_evals_FE(i) = evals_FE;
        h_list(i) = h_avg; %store actual average h

        %integration with explicit midpoint
        [~, X_list_EM, ~, evals_EM] = fixed_step_integration(rate_func_in, @explicit_midpoint_step, t_span, X0, h_ref);
        X_final_EM = X_list_EM(end, :)';
        error_EM(i) = norm(X_final_EM - X_analytical_final);
        num_evals_EM(i) = evals_EM;
    end

    %fit lines
    [p_FE_h, k_FE_h] = loglog_fit(h_list, error_FE);
    [p_EM_h, k_EM_h] = loglog_fit(h_list, error_EM);
    [p_FE_eval, k_FE_eval] = loglog_fit(num_evals_FE, error_FE);
    [p_EM_eval, k_EM_eval] = loglog_fit(num_evals_EM, error_EM);

    figure;
    loglog(h_list, error_FE, 'bo', 'DisplayName', 'Forward Euler Data');
    hold on;
    loglog(h_list, k_FE_h * h_list.^p_FE_h, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (p=%.2f)', p_FE_h));
    loglog(h_list, error_EM, 'ro', 'DisplayName', 'Explicit Midpoint Data');
    loglog(h_list, k_EM_h * h_list.^p_EM_h, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('EM Fit (p=%.2f)', p_EM_h));
    hold off;
    grid on;
    title(['Global Error vs. Step Size (Test Function ' num2str(flag) ')']);
    xlabel('Average Step Size (h)');
    ylabel('Global Error');
    legend('show', 'Location', 'southeast');

    figure;
    loglog(num_evals_FE, error_FE, 'bo', 'DisplayName', 'Forward Euler Data');
    hold on;
    loglog(num_evals_FE, k_FE_eval * num_evals_FE.^p_FE_eval, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (p=%.2f)', p_FE_eval));
    loglog(num_evals_EM, error_EM, 'ro', 'DisplayName', 'Explicit Midpoint Data');
    loglog(num_evals_EM, k_EM_eval * num_evals_EM.^p_EM_eval, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('EM Fit (p=%.2f)', p_EM_eval));
    hold off;
    grid on;
    title(['Global Error vs. Function Evals (Test Function ' num2str(flag) ')']);
    xlabel('Number of Rate Function Calls');
    ylabel('Global Error');
    legend('show', 'Location', 'southwest');
end

