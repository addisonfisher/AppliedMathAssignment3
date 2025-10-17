%Generalized local truncation solver
%Inputs: 
%rate_func_in: the function used to compute dxdt, to be used for the
%function used to calculate X+1
%t: the reference time used to calculate error across different solving
%methods 
%flag: denotes whether explicit_midpoint (1) or forward_euler_step (0) is
%used to calculate the error

function [p_FE, p_EM, p_analytical] = local_truncation(rate_func_in, t, flag)
    
    tref = t;
    h = logspace(-5, 1, 100);

    if flag == 1
        solution_func = @solution01;
    elseif flag == 2
        solution_func = @solution02;
    end

    %allocate arrays for errors
    error_FE = zeros(size(h));
    error_EM = zeros(size(h));
    analytical_diff = zeros(size(h));

    X_tref = solution_func(tref);
    
    for i = 1:length(h)
        h_i = h(i);
        
        %exact solution at t_ref + h
        X_tref_plus_h = solution_func(tref + h_i);
        
        %forward ruler
        G_x_FE = forward_euler_step(rate_func_in, tref, X_tref, h_i);
        error_FE(i) = norm(G_x_FE - X_tref_plus_h);

        %explicit midpoint
        G_x_EM = explicit_midpoint_step(rate_func_in, tref, X_tref, h_i);
        error_EM(i) = norm(G_x_EM - X_tref_plus_h);
        
        %|X(t+h) - X(t)| for comparison
        analytical_diff(i) = norm(X_tref_plus_h - X_tref);
    end

    %using loglog_fit to find p and k
    [p_FE, k_FE] = loglog_fit(h, error_FE);
    [p_EM, k_EM] = loglog_fit(h, error_EM);
    [p_analytical, k_analytical] = loglog_fit(h, analytical_diff);

    figure;
    set(gca, 'XScale', 'log', 'YScale', 'log')
    hold on;
    
    %plot collected data points
    loglog(h, error_FE, 'bo', 'DisplayName', 'Forward Euler Data');
    loglog(h, error_EM, 'ro', 'DisplayName', 'Explicit Midpoint Data');
    loglog(h, analytical_diff, 'ko', 'DisplayName', '|X(t+h)-X(t)| Data');

    %plot corresponding fit lines
    loglog(h, k_FE * h.^p_FE, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (p=%.2f)', p_FE));
    loglog(h, k_EM * h.^p_EM, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('EM Fit (p=%.2f)', p_EM));
    loglog(h, k_analytical * h.^p_analytical, 'k-', 'LineWidth', 2, 'DisplayName', sprintf('Analytical Diff Fit (p=%.2f)', p_analytical));
    
    hold off;
    title(['Local Truncation Error Analysis (Test Function ' num2str(flag) ')']);
    xlabel('Step Size (h)');
    ylabel('Error');
    legend('show', 'Location', 'northwest');
end