%Runs fixed step numerical integration
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%step_func: the function used to approximate X(t) at the next step
% [XB,num_evals] = step_func(rate_func_in,t,XA,h)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = ...
fixed_step_integration(rate_func_in,step_func,tspan,X0,h_ref)
    t_start = tspan(1);
    t_end = tspan(2);

    %get num steps over time interval
    num_steps = ceil((t_end - t_start) / h_ref);
    h_avg = (t_end - t_start) / num_steps;

    t_list = zeros(num_steps + 1, 1);
    %initialize x_list using number of ics in X0
    X_list = zeros(num_steps + 1, size(X0, 1)); 

    t_list(1) = t_start;
    X_list(1, :) = X0'; 
    
    total_num_evals = 0;

    for i = 1:num_steps
        current_t = t_list(i);
        %get current state vector
        current_X = X_list(i, :)'; 

        %get next state
        [next_X, step_evals] = step_func(rate_func_in, current_t, current_X, h_avg);
        
        %update the total count of rate function evaluations
        total_num_evals = total_num_evals + step_evals;

        t_list(i + 1) = t_list(i) + h_avg;
        X_list(i + 1, :) = next_X';
    end
    
    num_evals = total_num_evals;
end