



function [del1_ut,del2_ut,u,iterationTime] = directSolver_CollocationScheme_2D(N, k, m, R, theta, accuracy, q11_, q22_, epsilon, B1, B2, C1, C2, incField, j)

   tic        
   h1 = 2*pi/N; 
   h2 = (2*m*R+2*epsilon)/N;   
   alpha = k*cos(theta);  
   const = sqrt(2*pi*(2*m*R+2*epsilon));
   
% generate the grids
   [XN,YN] = meshgrid(-N/2+1:N/2);  %XN-left-right, YN-down-up        
        
% Fourier coefficients of the the smoothed kernel      
   KernelHat = kernel_hat2(R,m,k,alpha,XN,YN,epsilon);
  %KernelHat = cutoffKernel_hat2(R,m,k,N,alpha,epsilon);      
  
%   uij = inline('exp(1i*((alpha+j)*x1+sqrt(k^2-(alpha+j)^2)*x2)) + exp(1i*((alpha+j)*x1-sqrt(k^2-(alpha+j)^2)*x2))','x1','x2','k','alpha','j'); % incident field        
%   uij_ = uij(XN*h1,YN*h2,k,alpha,j);  % grid values of incident field  
 
    uij_down = exp(1i*((alpha + j)*XN*h1 - sqrt(k^2 - (alpha + j)^2)*YN*h2));    
    uij_up = exp(1i*((alpha + j)*XN*h1 + sqrt(k^2 - (alpha + j)^2)*YN*h2)); 
 
 if incField == 1          
     del1_incField = 1i*(alpha+j).*(uij_down + uij_up);  
     del2_incField = 1i*sqrt(k^2-(alpha+j)^2).*(uij_up - uij_down);     
 else       
     del1_incField = 1i*(alpha+j).*(uij_down - uij_up);  
     del2_incField = -1i*sqrt(k^2-(alpha+j)^2).*(uij_up + uij_down);      
 end

       
% we compute the right hand side
  RHS = (1i*(XN+alpha).*vfft2Alpha2(q11_.*del1_incField,B1,B2,N,alpha,const) +...
         1i*pi/(m*R+epsilon)*YN.*vfft2Alpha2(q22_.*del2_incField,B1,B2,N,alpha,const)).*const.*KernelHat;  
    
% GMRES iteration solver        
  u_hat = gmres(@(x)matrixHCollocation2(x,q11_,q22_,R,m,N,alpha,XN,YN,KernelHat,B1,B2,C1,C2,epsilon), RHS(:),100,accuracy,300);     
  iterationTime = toc;       
            
% compute solution from its Fourier coefficients
  u = ivfft2Alpha2(u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);    
  del1_ut = ivfft2Alpha2(1i*(XN+alpha).*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const) + del1_incField;  
  del2_ut = ivfft2Alpha2(1i*pi/(m*R+epsilon).*YN.*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const) + del2_incField;
   
   
   