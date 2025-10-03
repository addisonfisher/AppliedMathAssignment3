%This function computes the value of X at the next time step
%using the explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    num_evals = 0;
    
    %get slope from first "half-step" and apply eulers
    slope_initial = rate_func_in(t,XA);
    num_evals = num_evals + 1;

    X_mid = XA + (h / 2) * slope_initial;
    t_mid = t + h / 2;
    
    %get slope for second half and apply
    slope_midpoint = rate_func_in(t_mid, X_mid);
    XB = XA + h * slope_midpoint;
end