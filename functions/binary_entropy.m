function Hb=binary_entropy(p)
% Binary entropy function
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
    Hb=-p*log2(p)-(1-p)*log2(1-p);
return