function Pr_s=Get_MBPMF(s,rho)
% function Pr_s=Get_MBPMF(s,rho)
% generate Maxwell-Boltzmann distribution
%
%   Input:
%   - s(D,M) :     symbols of the M-ary modulation alphabet (M D-dimensional symbols)
%   - rho    :     coefficient of M-B distribution
%
%   Output:
%   - Pr_s(M):     probability mass function of constellation s
%   
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
tmp = exp(-rho*sum(s.^2,2));
Pr_s = (tmp./sum(tmp)).';
return