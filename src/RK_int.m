function int_value = RK_int(obj_func_mat, upper_limit, lower_limit, interval)
%%
% Runge-Kutta integral numerical method to obtain the integral values.
% input: 
%       obj_func_mat -- the integrated object function's values,
%                       corresponding to the range from upper to lower 
%                       integral limits;
%       upper_limit  -- the upper integral limit;
%       lower_limit  -- the lower integral limit;
%       interval     -- the integal interval between every two adjacent 
%                       object values in matrix "obj_func_mat".
% output:
%       int_value    -- the result of numerical integral.
%%



trapz()
