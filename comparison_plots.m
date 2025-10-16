function comparison_plots
    t_span = [0,10];
    X0 = 1;
    h_span = [0.1, 0.5];
    h_step = 0.1;

    h_ref = h_span(1):h_step:h_span(2);

    %Initializing values
    T1_FE = cell(1, length(h_ref));
    X1_BE = cell(1, length(h_ref));

    T1_BE = cell(1, length(h_ref));
    X1_BE = cell(1, length(h_ref));

    T1_EM = cell(1, length(h_ref));
    X1_EM = cell(1, length(h_ref));
    
    T1_IM = cell(1, length(h_ref));
    X1_IM = cell(1, length(h_ref));


    %Computing & Plotting the Forward Euler Approximation for various timesteps
    figure()
    for i =1:length(h_ref)
        [T1_FE{i}, X1_FE{i}, ~, ~] = fixed_step_integration(@rate_func01,@forward_euler_step,t_span,X0,h_ref(i));
        plot(T1_FE{i}, X1_FE{i},'--', 'DisplayName', 'Step Value ' + string(h_ref(i)));
        hold on;
    end
    X1_sol = solution01(T1_FE{1});
    plot(T1_FE{1}, X1_sol, 'r-', 'DisplayName', 'Solution');

    ylim([-1.5 1.5])
    legend();
    title('Forward Euler Approximations');
    hold off;

    %Computing & Plotting the Backward Euler Approximation for various timesteps
    figure()
    for i =1:length(h_ref)
        [T1_BE{i},  X1_BE{i}, ~, ~] = fixed_step_integration(@rate_func01, @backward_euler_step, t_span, X0, h_ref(i));
        plot(T1_BE{i}, X1_BE{i},'--', 'DisplayName', 'Step Value ' + string(h_ref(i)));
        hold on;
    end
    X1_sol = solution01(T1_BE{1});
    plot(T1_BE{1}, X1_sol, 'r-', 'DisplayName', 'Solution');
    ylim([-1.5 1.5]);
    legend();
    title('Backward Euler Approximations');
    hold off;

     %Computing & Plotting the Explicit Midpoint Approximation for various timesteps
    figure()
    for i =1:length(h_ref)
        [T1_EM{i}, X1_EM{i}, ~, ~] = fixed_step_integration(@rate_func01, @explicit_midpoint_step,t_span,X0,h_ref(i));
        plot(T1_EM{i}, X1_EM{i},'--', 'DisplayName', 'Step Value ' + string(h_ref(i)));
        hold on;
    end
    X1_sol = solution01(T1_EM{1});
    plot(T1_EM{1}, X1_sol, 'r-', 'DisplayName', 'Solution');

    ylim([-1.5 1.5])
    legend();
    title('Explicit Midpoint Approximations');
    hold off;

    %Computing & Plotting the Implicit Midpoint Approximation for various timesteps
    % figure()
    % for i =1:length(h_ref)
    %     [T1_IM{i}, X1_IM{i}, ~, ~] = fixed_step_integration(@rate_func01, @implicit_midpoint_step,t_span,X0,h_ref(i));
    %     plot(T1_IM{i}, X1_IM{i},'--', 'DisplayName', 'Step Value ' + string(h_ref(i)));
    %     hold on;
    % end
    % X1_sol = solution01(T1_IM{1});
    % plot(T1_IM{1}, X1_sol, 'r-', 'DisplayName', 'Solution');
    % 
    % ylim([-1.5 1.5])
    % legend();
    % title('Implicit Midpoint Approximations');
    % hold off;

end
