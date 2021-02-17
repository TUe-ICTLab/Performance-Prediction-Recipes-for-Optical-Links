function y = channel_phase_noise(x,sz2,sp2)
%function y = channel_phase_noise(x,sz2,sp2)
%   Simulate a channel with AWGN and memoryless phase noise
%   y and x are real 2xN matrices (2D-dimensional symbols, corresponding to complex symbols)
%
%   The complex AWGN noise has a total variance sz2.
%   The phase noise has a variance sp2
%
% E. Agrell, M. Secondini, A. Alvarado and T. Yoshida
% Feb. 2021

[D,N]=size(x);
if (D~=2)
    error('Error! The channel takes only 2D input symbols');
end

%% Convert 2D real symbols into complex symbols
xc=x(1,:)+1i*x(2,:);
%% Channel (complex input and complex output)
z=sqrt(sz2/2)*(randn(1,N)+1i*randn(1,N)); % AWGN noise
theta=sqrt(sp2)*randn(1,N);             % Phase noise
yc=xc.*exp(-1i*theta)+z;
%% Go back to 2D real vectors
y=[real(yc);imag(yc)];

end

