%Computes the theta derivatives of each vertex coordinate for the Jansen linkage
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% these are assumed to be legal values that are roots of the error funcs!
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%dVdtheta: a column vector containing the theta derivates of each vertex coord
function dVdtheta = compute_velocities(vertex_coords, leg_params, theta)

    link_error_func = @(coords) link_length_error_func(coords, leg_params);

    J = approximate_jacobian(link_error_func, vertex_coords);
    
    %compute d1 velocities (determined by crank angle)
    dx1_dtheta = -leg_params.crank_length * sin(theta);
    dy1_dtheta =  leg_params.crank_length * cos(theta);

    dx2_dtheta = 0;
    dy2_dtheta = 0;

    %total number of variables is 2 * num_vertices
    num_vars = 2 * leg_params.num_vertices;
    
    %m_top = first equations for v1 and v2
    M_top = [eye(4), zeros(4, num_vars - 4)];
    M = [M_top; J];
     
    %constructing b vector
    B = [dx1_dtheta; dy1_dtheta; dx2_dtheta; dy2_dtheta; zeros(leg_params.num_linkages, 1)];

    dVdtheta = M \ B;
end