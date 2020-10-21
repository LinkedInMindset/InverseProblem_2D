

 
% this function is to compute "vainikko inverse discrete fourier transform "in 2D
function y = ivfft22(x,C1,C2,N,const)
%const = sqrt(2*pi*(2*m*R+2*epsilon));
y = N^2./const.*ivfft1(ivfft1(x,C1,C2).',C1,C2).';
end

