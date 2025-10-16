% Example: simple exponential decay dx/dt = -2x
function comparison_plots
    t_span = [0,10];
    X0 = 1;
    h_ref = 0.1;

    [t1_approx1, X1_approx1, h1_avg_1, n1_evals_1] = forward_euler_fixed_step_integration(@rate_func01, t_span, X0, h_ref);
    [t1_approx2, X1_approx2, h1_avg_2, n1_evals_2] = explicit_midpoint_fixed_step_integration(@rate_func01, t_span, X0, h_ref);
    X1_sol = solution01(t1_approx1);

    figure();
    plot(t1_approx1, X1_sol, 'r-');
    hold on;
    plot(t1_approx1, X1_approx1, 'b--');

    legend('Solution', 'Forward Euler Approx');
    title('Forward Euler Integration');

    figure();
    plot(t1_approx2, X1_sol, 'r-');
    hold on;
    plot(t1_approx2, X1_approx2, 'b--');
    ylim([-1.5 1.5])
    legend('Solution', 'Explicit Midpoint Approx');

    xlabel('t'); ylabel('X(t)');
    title('Explicit Midpoint Integration');
end
