function day_work()
    X0 = [0.8308;0.5853;0.5497];%randn(3,1);

    [~, J_analytical] = test_function01(X0);
    J_numerical = approx_jacobian01(@test_function01,X0);

    solver_params = struct();
    solver_params.dxmin = 1e-10;
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.maxiter = 200;
    solver_params.approx = 1;

    Xguess = rand(3,1);

    xroot = multinewt(@test_function01,Xguess,solver_params);

    disp(xroot)

    f_root = test_function01(xroot);

    disp(f_root)

    % disp(J_analytical)
    % disp(J_numerical)

end


function xroot = multinewt(fun,X,solverparams)

    dxmin = solverparams.dxmin;
    ftol = solverparams.ftol;
    dxmax = solverparams.dxmax;
    maxiter = solverparams.maxiter;
    approx = solverparams.approx;


    if approx
        fval = fun(X);
        J = approx_jacobian01(fun,X)
    else
        [fval,J] = fun(X);
    end

    delta_x = 1;

    count = 0;

    while count<maxiter && norm(delta_x)>dxmin && norm(fval) > ftol && norm(delta_x)<dxmax
        count = count+1;

        if approx
        fval = fun(X);
        J = approx_jacobian01(fun,X);
        else
        [fval,J] = fun(X);
        end

        delta_x = -J\fval;

        X = X + delta_x;





    end
    xroot = X
end


function J = approx_jacobian01(fun,X)

    f0 = fun(X);
    J = zeros(length(f0), length(X));


    e_n = zeros(length(X),1);

    delta_X = 1e-6;

    for n = 1:length(X)
        e_n(n) = 1;


        f_plus = fun(X+e_n*delta_X);
        f_minus = fun(X-e_n*delta_X);

        J(:,n) = (f_plus-f_minus) / (2*delta_X);

        e_n(n) = 0;
    end

end

function [fval, J] = test_function01(X)
    x1 = X(1); x2 = X(2); x3 = X(3);

    
    f1 = x1^2 + x2^2 - 6 + x3^5;
    f2 = x1*x3 + x2 -12;
    f3 = sin(x1 + x2 + x3);


    

    fval = [f1;f2;f3];

    J = [ [2*x1, 2*x2, -5*x3^4];...
        [x3, 1, x1];...
        [cos(x1+x2+x3),cos(x1+x2+x3),cos(x1+x2+x3)]];
end