function test()  
    clear;
    clc;
    close all;
    
    leg_params = generate_params();
    
    vertex_coords_test = [
         0;  50; % Vertex 1 guess
       -50;   0; % Vertex 2 guess
       -50;  50; % Vertex 3 guess
      -100;   0; % Vertex 4 guess
      -100; -50; % Vertex 5 guess
       -50; -50; % Vertex 6 guess
       -50;-100  % Vertex 7 guess
    ];
    
    figure;
    hold on;
    xlabel('x position');
    ylabel('y position');
    
    leg_drawing = initialize_leg_drawing(leg_params);
    
    update_leg_drawing(vertex_coords_test, leg_drawing, leg_params);
    
    axis equal;
    grid on;
end