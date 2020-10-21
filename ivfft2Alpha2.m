



% this function is to compute "alpha inverse discrete fourier transform "in 2D
function y = ivfft2Alpha2(x,C1,C2,N,alpha,const)
C = exp(1i.*alpha.*meshgrid(-N/2+1:N/2).*2.*pi/N);
y = C.*ivfft22(x,C1,C2,N,const); 
end