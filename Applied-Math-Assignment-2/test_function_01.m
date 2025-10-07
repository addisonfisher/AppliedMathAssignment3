%Inputs values of X

function [f_val, J] = test_function_01(X)
    x1 = X(1); x2 = X(2); x3 = X(3);

    
    f1 = x1^2 + x2^2 - 6 + x3^5;
    f2 = x1*x3 + x2 -12;
    f3 = sin(x1 + x2 + x3);


    

    f_val = [f1;f2;f3];

    J = [ [2*x1, 2*x2, -5*x3^4];...
        [x3, 1, x1];...
        [cos(x1+x2+x3),cos(x1+x2+x3),cos(x1+x2+x3)]];
end