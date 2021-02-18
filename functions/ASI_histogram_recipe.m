function ASItmp = ASI_histogram_recipe(aLsq,Bin)
% function ASItmp = ASI_histogram_recipe(aLsq,Bin)
% ASI computation based on normalized and quantized asymmetric L-value
%
%   Input:
%   - aLsq   :    assymmetric (scaled and quantized) L-value
%   - Bin    :    number of bins in L-value quantization
%   
%   Output:
%   - ASItmp :    output asymmetric information
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021
halfBin = (Bin-mod(Bin,2))/2;
Lambda_pos = zeros(halfBin,1);
Lambda_neg = zeros(halfBin,1);
sum_tmp = zeros(halfBin,1);
for x=1:1:halfBin
    Lambda_pos(x,1) = length(aLsq(aLsq == -x))/length(aLsq);
    Lambda_neg(x,1) = length(aLsq(aLsq == x))/length(aLsq);
    if Lambda_pos(x,1) ~= 0
        sum_tmp(x,1) = Lambda_pos(x,1)*log2(2*(Lambda_pos(x,1))/(Lambda_pos(x,1)+Lambda_neg(x,1)));
    end
    if Lambda_neg(x,1) ~= 0
        sum_tmp(x,1) = sum_tmp(x,1) + Lambda_neg(x,1)*log2(2*Lambda_neg(x,1)/(Lambda_pos(x,1)+Lambda_neg(x,1)));
    end
end
ASItmp = sum(sum_tmp);
return