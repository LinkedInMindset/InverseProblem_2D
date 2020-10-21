

% This function computes Fourier coefficients of the smoothed quasi-periodic kernel


function y = cutoffKernel_hat2(R,m,k,N,alpha,epsilon)


% R =rho

% n = 8; % degree of the approximating polynomial

[X,Y] = meshgrid(-N/2+1:N/2); 

%[X2N,Y2N] = meshgrid(-2*N/2+1:2*N/2);

y = zeros(N); 

p1 = polynomial_hat2(2*N,m,R,epsilon);

p2 = kernel_hat2(R,m,k,alpha,X,Y,epsilon);

const = sqrt(2*pi*(2*m*R+2*epsilon));


%for j1 = -1*N/2+1:1*N/2
    
   %for  j2 = -1*N/2+1:1*N/2

    for l2 = -N/2+1:N/2                   
                                
        % y(j2+N/2,j1+N/2) = y(j2+N/2,j1+N/2) + p2(l2 + N/2 + N*(j1 + N/2 - 1)).*p1(j2 - l2 + N + 2*N*(j1 + N - 1))/sqrt(8*pi*R);
        
        % y(j2+N/2,j1+N/2) = y(j2+N/2,j1+N/2) + p2(l2 + N/2 ,j1 + N/2).*p1(j2 - l2 + N ,0 + N)/const;
         
        %  y(j2+N/2,j1+N/2) =  sum(p2((-N/2+1:N/2) + N/2 ,j1 + N/2).*p1(j2 - (-N/2+1:N/2) + N ,0 + N)/const);
              
        % y = y + p2(l2 + 2*N/2 + 2*N*(X + 2*N/2 - 1)).*p1(Y - l2 + 3*N/2 + 3*N*(0 + 3*N/2 - 1))/sqrt(2*pi*(2*m*R+2*epsilon));   
         
         y = y + p2(l2 + N/2 + N*(X + N/2 - 1)).*p1(Y - l2 + N + 2*N*(0 + N - 1))/const;
         
    end    
    
   %end
   
%end

 
%  mesh((-N/2+1:N/2),(-N/2+1:N/2),real(y));
%  axis square
%  colorbar
