function stability_comparison_plots()
    % set up parameters
    t_span = [0, 20];
    x0 = 1;
    h_refs = [0.38, 0.45];
    rate_func = @rate_func01;
    analytical_sol = @solution01;

    % method handles and names
    step_funcs = {@forward_euler_step, @backward_euler_step, @explicit_midpoint_step, @implicit_midpoint_step};
    method_names = {'Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint'};

    % loop through h_ref values
    for h = h_refs
        figure;
        sgtitle(['Comparison for h_{ref} = ' num2str(h)]);

        % loop through methods
        for i = 1:length(step_funcs)
            % compute numerical solution
            [t_num, x_num, ~, ~] = fixed_step_integration(rate_func, step_funcs{i}, t_span, x0, h);

            % compute analytical solution
            x_analytical = analytical_sol(t_num);

            subplot(4, 1, i);
            plot(t_num, x_num, 'b-', 'DisplayName', 'Numerical');
            hold on;
            plot(t_num, x_analytical, 'r--', 'DisplayName', 'Analytical');
            hold off;

            % add labels and title
            title(method_names{i});
            xlabel('time (t)');
            ylabel('x(t)');
            legend;
            grid on;
        end
    end
end