%This function computes the value of X at the next time step
%using the Backward Euler approximation
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
function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
    %define func to get root by passing through newton's
    G = @(X_next_guess) XA + h * rate_func_in(t + h, X_next_guess) - X_next_guess;
    
    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.maxiter = 200;
    solver_params.approx = 1;

    X_initial = XA;

    [XB, num_evals, ~] = multi_newton_solver(G, X_initial, solver_params);

end