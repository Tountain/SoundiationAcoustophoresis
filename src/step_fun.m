function y = step_fun(T)
%%
% input a matrix 'T', and return a matrix the same length as 'T' with zeros
% where 'abs(T) > 1' and ones where 'abs(T) <= 1'.
% Ref. stepfun.m of the incode function in MATLAB.
%%

[n, m] = size(T);
y = ones(n, m);
[rr, cc, ~] = find(T > 1);
for ii = 1:length(rr)
    y(rr(ii), cc(ii)) = 0;
end

%%