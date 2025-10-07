%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function [J, num_evals] = approximate_jacobian(fun,x)
    num_inputs = length(x); 
    %evaluate fun once to get the number of outputs
    f_val = fun(x);
    num_outputs = length(f_val);
    
    %initialize jacobian
    J = zeros(num_outputs, num_inputs);

    %setting step size
    h = 1e-6;

    for i = 1:num_inputs
        %creating e matrices-- [1;0...0], [0;1;0...0]
        e = zeros(num_inputs, 1);
        e(i) = 1;

        J(:, i) = (fun(x + h*e) - fun(x - h*e)) / (2*h);
    end

    num_evals = 2 * num_inputs;

end