function Hk_n_p = sHankel_partial(n,x,type)
%%
% Obtain the differentiate of Spherical Hankel Function
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
% Based on Eq.(6.69) at E.G.W@1999
%%

Hk_n_p = sHankel(n-1,x,type) - (n+1)*sHankel(n,x,type)/x;

