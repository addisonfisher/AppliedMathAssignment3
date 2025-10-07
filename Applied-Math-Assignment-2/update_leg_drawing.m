% Updates the plot objects that visualize the leg linkage
% for the current leg configuration
% INPUTS:
% complete_vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params)
    
    % Convert the column vector of coordinates into a more usable n-by-2 matrix
    vertex_matrix = column_to_matrix(complete_vertex_coords);

    for linkage_index = 1:leg_params.num_linkages
        % Get the indices of the two vertices connected by the current link
        vertex1_index = leg_params.link_to_vertex_list(linkage_index, 1);
        vertex2_index = leg_params.link_to_vertex_list(linkage_index, 2);

        % Get the coordinates of these two vertices from the matrix
        x1 = vertex_matrix(vertex1_index, 1);
        y1 = vertex_matrix(vertex1_index, 2);
        x2 = vertex_matrix(vertex2_index, 1);
        y2 = vertex_matrix(vertex2_index, 2);
        
        % Create two-element arrays for the x and y data of the line
        line_x = [x1, x2];
        line_y = [y1, y2];
        
        set(leg_drawing.linkages{linkage_index}, 'xdata', line_x, 'ydata', line_y);
    end
    for vertex_index = 1:leg_params.num_vertices
        dot_x = vertex_matrix(vertex_index, 1);
        dot_y = vertex_matrix(vertex_index, 2);
        
        set(leg_drawing.vertices{vertex_index}, 'xdata', dot_x, 'ydata', dot_y);
    end

    crank_origin = leg_params.vertex_pos0;
    crank_tip_coords = vertex_matrix(1, :);
    
    crank_x = [crank_origin(1), crank_tip_coords(1)];
    crank_y = [crank_origin(2), crank_tip_coords(2)];

    set(leg_drawing.crank, 'xdata', crank_x, 'ydata', crank_y);
end