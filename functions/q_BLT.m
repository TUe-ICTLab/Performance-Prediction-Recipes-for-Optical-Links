function q=q_BLT(y,x,sz2,sp2)
% Function q given using Eq. (3) of 
% T. Koike-Akino, D. S. Millar, K. Kojima, and K. Parsons, “Phase noise-
% robust LLR calculation with linear/bilinear transform for LDPC-coded
% coherent communications,” in Proc. Conference on Lasers and Electro-
% Optics (CLEO), San Jose, CA, May 2015. 
%
% For more details, see footnote 2 in 
% "Performance Prediction Recipes for Optical Links", Photonics Technology
% Letters, 2021, by Agrell, Secondini, Alvarado and Yoshida.
%
% In this function y and x are real 2xN matrices (2D-dimensional symbols,
% to be converted into complex symbols) 
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

yc=y(1,:)+1i*y(2,:);
xc=x(1,:)+1i*x(2,:);
R=(3*yc-xc)/2;
S=(yc+xc)/2;
q=exp(-abs(R-S).^2/sz2+(2*sp2/sz2)./(sz2+2*sp2*abs(S).^2).*(imag(conj(S).*R)).^2-1/2*log(sz2+2*sp2*abs(S).^2));   
return