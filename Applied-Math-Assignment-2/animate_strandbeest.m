function animate_strandbeest()

    close all;
    leg_params = generate_params();
    

    global dxtip_analytical;
    global dytip_analytical;
    
    figure;
    hold on;
    axis equal;
    axis([-120 40 -120 60]); 
    title('Strandbeest Linkage Animation');
    xlabel('x position');
    ylabel('y position');
    grid on;
    
    leg_drawing = initialize_leg_drawing(leg_params);
    vertex_coords_guess = [
        0; 50; ...      % Vertex 1
        -50; 0; ...     % Vertex 2
        -50; 50; ...    % Vertex 3
        -100; 0; ...    % Vertex 4
        -100; -50; ...  % Vertex 5
        -50; -50; ...   % Vertex 6
        -50; -100      % Vertex 7
    ];

    leg_tip_history = [];

    num_rotations = 10;
    theta_steps = linspace(0, num_rotations * 2 * pi, (num_rotations*200)); 

    trace_line = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);

    count = 1;

    q = [];
    
    for theta = theta_steps
        vertex_coords_root = compute_coords(vertex_coords_guess, leg_params, theta);
        
        update_leg_drawing(vertex_coords_root, leg_drawing, leg_params);

        leg_tip_x = vertex_coords_root(13);
        leg_tip_y = vertex_coords_root(14);
        
        
        q = quiver(leg_tip_x,leg_tip_y,dxtip_analytical(count),dytip_analytical(count),'r','LineWidth',1,'AutoScale','on');
        count = count+1;

        if count >= 200
            count = 1;
        end

        leg_tip_history(:, end+1) = [leg_tip_x; leg_tip_y];
        set(trace_line, 'xdata', leg_tip_history(1,:), 'ydata', leg_tip_history(2,:));
        
        vertex_coords_guess = vertex_coords_root;
        
        drawnow;
        pause(0.01);

        delete(q);
    end
end