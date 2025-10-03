function midpoint_verification()
    rate_func = @rate_func01;      
    tspan = [0, 15];               
    X0 = 1;                        
    h_refs = [0.5, 0.2, 0.05];     
    
    figure;
    hold on;
    grid on;
    title('Explicit Midpoint Verification');
    xlabel('Time (t)');
    ylabel('State (x)');
    
    legend_entries = {};
    for i = 1:length(h_refs)
        h = h_refs(i);
        [t_num, X_num, ~, ~] = explicit_midpoint_fixed_step_integration(rate_func, tspan, X0, h);
        plot(t_num, X_num, 'LineWidth', 1.5);
        legend_entries{end+1} = ['Numerical, h_{ref} = ' num2str(h)];
    end
    
    t_analytical = linspace(tspan(1), tspan(2), 500);
    X_analytical = solution01(t_analytical); 
    
    plot(t_analytical, X_analytical, 'k--', 'LineWidth', 2);
    legend_entries{end+1} = 'Analytical Solution (cos(t))';
    
    legend(legend_entries, 'Location', 'best');
    hold off;
end