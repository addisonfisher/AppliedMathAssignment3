%Runs numerical integration using forward Euler approximation

%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)

%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration

function [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
    
    t_list = [];
    X_list = [];
    X_list(end) = X0;
    h_list = [];
    t_start = t_span(1);
    eval_num = [];

    %Looping through t + hvals while less than t_span(2)
        h_list(end) = t_val - t_start;

        [X_list(end), eval_num(end)] = forward_euler_step(rate_func_in,t_val,X_list(end),h);

        t_start = T_val;

    %end of loop
    
    h_avg = mean(h_list);
    num_evals = sum(h_avg);

end