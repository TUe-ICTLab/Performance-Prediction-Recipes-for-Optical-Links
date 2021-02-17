function q=q_awgn(y,x,sigma2)
% Function q in (8)
% y and x are real DxN matrices (N D-dimensional symbols)
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

q=exp(-sum(abs(y-x).^2,1)/(2*sigma2));

return