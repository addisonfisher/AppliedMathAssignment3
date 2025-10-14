%Generalized local truncation solver
%Inputs: 
%rate_func_in: the function used to compute dxdt, to be used for the
%function used to calculate X+1
%t: the reference time used to calculate error across different solving
%methods 
%flag: denotes whether explicit_midpoint (1) or forward_euler_step (0) is
%used to calculate the error

function p = local_truncation(rate_func_in, t, flag)
    
    tref = t;
    h = logspace(-5, 1, 100);
    error_list = [];

    for i = 1: length(h)

        if flag == 1
            X = solution01(tref);
            X_h = solution01(tref+h(i)); %make sure that h is a singular value when debugging
            G_x = explicit_midpoint_step(rate_func_in, tref, X, h);
        elseif flag == 0
            X = solution02(tref);
            X_h = solution02(tref+h(i));
            G_x = forward_euler_step(rate_func_in, tref, X, h);
        end 

        error = norm(G_x - X_h);
        error_list(end + 1) = error; 
    end
    
    [p,~] = loglog_fit(h, error_list);

end
