function [q,ea,iter] = romberg(f,a,b,es,maxIt)
%romberg: Romber integration quadrature
% q = romberg(f,a,b,es,maxIt): Romberg Integration
%
% inputs:
% f = function to be integrated
% a,b = integration limits
% es = desired relative error(default = 0.0000001%)
% maxIt = maximum allowable iterations (default = 50)
% outputs:
% q = integral estimate
% ea = magnitude of percent approximate relative error
% iter = number of iterations

% Created by: David Pelley
% 9/26/2019

if nargin<3,error('At least 3 input arguments required'), end
if nargin<4||isempty(es), es = 0.000001; end
if nargin<5||isempty(maxIt), maxIt = 50; end

I(1,1) = trap(f,a,b,1);
iter = 0;
while iter<maxIt
    iter = iter+1;
    n = 2^iter;
    I(iter+1,1) = trap(f,a,b,n);
    for k = 2:iter+1
        j = 2+iter-k;
        I(j,k) = (4^(k-1)*I(j+1,k-1)-I(j,k-1))/(4^(k-1)-1);
    end
    ea = abs((I(1,iter+1)-I(2,iter))/I(1,iter+1))*100;
    if ea < es, break; end
end
q = I(1,iter+1);