function Bsl_p = sBessel_partial(n,x,type)
%%
% Obtain the differentiate of Spherical Bessel Function 
% n--order
% x--independent variable (SHOULD be row vector)
% type--the first kind or the second kind
% Based on Eq.(6.69) at E.G.W@1999
%%

Bsl_p = sBessel(n-1,x,type) - (n+1)*sBessel(n,x,type)/x;

