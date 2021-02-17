function Hb=binary_entropy(p)
    % Binary entropy function
    Hb=-p*log2(p)-(1-p)*log2(1-p);
return