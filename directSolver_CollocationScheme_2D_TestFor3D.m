



function [u iterationTime] = directSolver_CollocationScheme_2D_TestFor3D(N, k, m, R, theta, q, accuracy, epsilon, j)

   tic        
   h1 = 2*pi/N; 
   h2 = (2*m*R+2*epsilon)/N;   
   alpha = k*cos(theta);  
   alphaj = alpha+j;
   betaj = sqrt(k^2-alphaj^2);
   const = sqrt(2*pi*(2*m*R+2*epsilon));
   
   
  B1=ones(N,1)*exp(-2*pi*1i*(1-N/2)*(1:N)/(N));
  B2=ones(N,1)*exp(-2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
  C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/(N));
  C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
   
% generate the grids
   [XN,YN] = meshgrid(-N/2+1:N/2);  % XN-left-right, YN-down-up     
   
   q_ = zeros(N);  
   q_(abs(YN*h2)<=R) = q;
        
% Fourier coefficients of the the smoothed kernel      
   KernelHat = kernel_hat2(R,m,k,alpha,XN,YN,epsilon);
  %KernelHat = cutoffKernel_hat2(R,m,k,N,alpha,epsilon);      
  
   %uij = inline('alpha/k*exp(1i*((alpha+j)*x1+sqrt(k^2-(alpha+j)^2)*x2)) - alpha/k*exp(1i*((alpha+j)*x1-sqrt(k^2-(alpha+j)^2)*x2))','x1','x2','k','alpha','j'); % incident field       
   %uij_ = uij(XN*h1,YN*h2,k,alpha,j);  % grid values of incident field  
   
   
   PlaneWave_down = alphaj/sqrt(abs(betaj)^2 + alphaj^2)*exp(1i*(alphaj*XN*h1 - betaj*YN*h2));  
     %PlaneWave_up = cos(theta)*exp(1i*(alphaj*XN*h1 + betaj*YN*h2));  
     %PlaneWave_up = zeros(N);
       
% we compute the right hand side
  RHS = (1i*(XN+alpha).*vfft2Alpha2(q_.*1i.*alphaj.*(PlaneWave_down),B1,B2,N,alpha,const) +...
         1i*pi/(m*R+epsilon)*YN.*vfft2Alpha2(q_.*1i.*(-betaj).*(PlaneWave_down),B1,B2,N,alpha,const)).*const.*KernelHat;  
    
% GMRES iteration solver        
  u_hat = gmres(@(x)matrixHCollocation2(x,q_,R,m,N,alpha,XN,YN,KernelHat,B1,B2,C1,C2,epsilon), RHS(:),30,accuracy,300);     
  iterationTime = toc;       
            
% compute solution from its Fourier coefficients
  u = ivfft2Alpha2(u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);    
 % del1_u = ivfft2Alpha2(1i*(XN+alpha).*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);  
 % del2_u = ivfft2Alpha2(1i*pi/(m*R+epsilon).*YN.*u_hat(YN + N/2 + N*(XN + N/2-1)),C1,C2,N,alpha,const);
   
 
 mesh((-N/2+1:N/2)*h1,(-N/2+1:N/2)*h2,real(u));  
 
    set(gca,'YDir','normal')    
    axis tight
    colorbar
    shading interp 