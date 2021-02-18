function Hs=entropy(Ps)
% function Hs=entropy(Ps)
% Entropy function
% 
%   Input:
%   - Ps(M)  :    probability mass function of transmitted symbol (M-ary over D-dimensions)
%   
%   Output:
%   - Hs     :    symbol entropy
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
Hs=sum(-Ps.*log2(Ps),2);
return
