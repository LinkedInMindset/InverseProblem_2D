



% this function is to compute A*x, where A is the left hand side matrix.

function x = matrixHCollocation2(x,q11_,q22_,R,m,N,alpha,XN,YN,khat,B1,B2,C1,C2,epsilon)  

const = sqrt(2*pi*(2*m*R+2*epsilon));    

x = x(YN+N/2+N*(XN+N/2-1));  % transfer vector x into matrix form 

x = x - (1i*(XN+alpha).*vfft2Alpha2(q11_.*ivfft2Alpha2(1i*(XN+alpha).*x,C1,C2,N,alpha,const),B1,B2,N,alpha,const) + ...
    1i*pi/(m*R+epsilon)*YN.*vfft2Alpha2(q22_.*ivfft2Alpha2(1i*pi/(m*R+epsilon)*YN.*x,C1,C2,N,alpha,const),B1,B2,N,alpha,const)).*const.*khat;

x = x(:);  % transfer the result (matrix form) into vector form

end