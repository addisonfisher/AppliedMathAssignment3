%Error function that encodes all necessary linkage constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%error_vec: a vector describing each constraint on the linkage
% when error_vec is all zeros, the constraints are satisfied
function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    distance_errors = link_length_error_func(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end