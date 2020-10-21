



function [del1_u del2_u u] = directSolver_GalerkinScheme_2D(N, k, m, R, theta, accuracy, q_3N, epsilon, B1, B2, C1, C2, X3N, Y3N, j)

tic

  h1 = 2*pi/N; 
  h2 = (2*m*R + 2*epsilon)/N;   
  alpha = k*cos(theta);  
  const = sqrt(2*pi*(2*m*R+2*epsilon));           
  
 % generate the grids
  [XN,YN] = meshgrid(-N/2+1:N/2);  % XN-left-right, YN-down-up        
   
% Fourier coefficients of the the smoothed kernel      
   KernelHat = kernel_hat2(R,m,k,alpha,XN,YN,epsilon);
  %KernelHat = cutoffKernel_hat2(R,m,k,N,alpha,epsilon);   
         
   
  uij = inline('exp(1i*((alpha+j)*x1+sqrt(k^2-(alpha+j)^2)*x2)) + exp(1i*((alpha+j)*x1-sqrt(k^2-(alpha+j)^2)*x2))','x1','x2','k','alpha','j'); % incident field
  
  uij_ = uij(X3N*h1/3,Y3N*h2/3,k,alpha,j);  % grid values of incident field 
 
  uij_down = exp(1i*((alpha + j)*X3N*h1/3 - sqrt(k^2-(alpha + j)^2)*Y3N*h2/3));
  
  uij_up = exp(1i*((alpha + j)*X3N*h1/3 + sqrt(k^2-(alpha + j)^2)*Y3N*h2/3));
       
% we compute the right hand side
 
  RHS = (1i.*(XN + alpha).*Res(vfft2Alpha2(q_3N.*uij_.*1i*(alpha + j),B1,B2,3*N,alpha,const),3*N,N) +...
         1i*pi/(m*R+epsilon).*YN.*Res(vfft2Alpha2(q_3N.*1i.*sqrt(k^2 - (alpha + j)^2).*(uij_up - uij_down),B1,B2,3*N,alpha,const),3*N,N)).*const.*KernelHat;  
    
% GMRES iteration solver        
   u_hat = gmres(@(x)matrixHGalerkin2(x,q_3N,R,m,N,alpha,XN,YN,KernelHat,B1,B2,C1,C2,epsilon), RHS(:),30,accuracy,300);  
   
%  iterationTime = toc;       
            
% compute solution from its Fourier coefficients
   C1 = ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/N);   
   C2 = ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+N^2/4-1))/N);
   
   
   u = ivfft2Alpha2(u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);    
   del1_u = ivfft2Alpha2(1i*(XN+alpha).*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);  
   del2_u = ivfft2Alpha2(1i*pi/(m*R+epsilon).*YN.*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);
   
   
   