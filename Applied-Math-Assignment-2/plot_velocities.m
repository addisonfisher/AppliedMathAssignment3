%Plots the velocities
function plot_velocities()
    leg_params = generate_params();
    
    theta_values = linspace(0, 2 * pi, 200);
    global dxtip_analytical;
    global dytip_analytical;
    dxtip_analytical = zeros(size(theta_values));
    dytip_analytical = zeros(size(theta_values));
    dxtip_numerical = zeros(size(theta_values));
    dytip_numerical = zeros(size(theta_values));
    
    vertex_coords_guess = [
        0; 50; ...
        -50; 0; ...
        -50; 50; ...
        -100; 0; ...
        -100; -50; ...
        -50; -50; ...
        -50; -100];

    for i = 1:length(theta_values)
        theta = theta_values(i);
        
        %vertex coordinates for the current theta
        vertex_coords_root = compute_coords(vertex_coords_guess, leg_params, theta);
        
        %analytical vel calc
        dVdtheta_analytical = compute_velocities(vertex_coords_root, leg_params, theta);
        dxtip_analytical(i) = dVdtheta_analytical(13);
        dytip_analytical(i) = dVdtheta_analytical(14);
        
        %numerical velocity calc
        compute_coords_func = @(t) compute_coords(vertex_coords_root, leg_params, t);
        J_numerical = approximate_jacobian(compute_coords_func, theta);
        dxtip_numerical(i) = J_numerical(13);
        dytip_numerical(i) = J_numerical(14);
        
        vertex_coords_guess = vertex_coords_root;
    end
    
    figure();
    plot(theta_values, dytip_analytical, 'b-', 'LineWidth', 2);
    hold on;
    plot(theta_values, dytip_numerical, 'r--', 'LineWidth', 2);
    title('Comparison of dy_{tip}/d\theta');
    xlabel('\theta (radians)');
    ylabel('dy_{tip}/d\theta');
    legend('Analytical', 'Numerical');
    grid on;

    figure();
    plot(theta_values, dxtip_analytical, 'b-', 'LineWidth', 2);
    hold on;
    plot(theta_values, dxtip_numerical, 'r--', 'LineWidth', 2);
    title('Comparison of dx_{tip}/d\theta');
    xlabel('\theta (radians)');
    ylabel('dx_{tip}/d\theta');
    legend('Analytical', 'Numerical');
    grid on;
end