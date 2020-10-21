



% this function is to compute A*x, where A is the left hand side matrix.
function x = matrixHGalerkin2(x,q_,R,m,N,alpha,X,Y,khat,B1,B2,C1,C2,epsilon)
  
%[X,Y]=meshgrid(-N/2+1:N/2); 
const = sqrt(2*pi*(2*m*R+2*epsilon));
x = x(Y+N/2+N*(X+N/2-1));  % transfer vector x into matrix form 
x = x - (1i.*(X+alpha).*Res(vfft2Alpha2(q_.*ivfft2Alpha2(Ext(1i.*(X+alpha).*x,N,3*N),C1,C2,3*N,alpha,const),B1,B2,3*N,alpha,const),3*N,N) + ...
    1i.*pi/(m*R+epsilon).*Y.*Res(vfft2Alpha2(q_.*ivfft2Alpha2(Ext(1i.*pi./(m*R + epsilon).*Y.*x,N,3*N),C1,C2,3*N,alpha,const),B1,B2,3*N,alpha,const),3*N,N)).*const.*khat;
x = x(:);  % transfer the result (matrix form) into vector form
end