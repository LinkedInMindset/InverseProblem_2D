


% this function is to compute "inverse vainikko discrete fourier transform "in 1D

function y=ivfft1(x,C1,C2)
%C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/N);
%u=ifft(x.*C1.');
%C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+N^2/4-1))/N);
y=ifft(x.*C1.').*C2.';
end
