function X=discreteinvrnd(p,m,n)
% function X=discreteinvrnd(p,m,n)
% generate variable following predetermined distribution
%
%   Input:
%   - p(M)   :    probability
%   - m      :    number of columns of the output matrix
%   - n      :    number of rows of the output matrix
%   
%   Output:
%   - X(m,n) :    output data
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
X = zeros(m,n);
for i = 1:m*n
    u = rand;
    I = find(u < cumsum(p));
    X(i) = min(I);
end
return