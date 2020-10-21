

% this function is to compute the near field using the direct data generating from "directSolver_GalerkinScheme_2D"

function   NearField = data_GalerkinScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q_N, q_3N, epsilon, B1, B2, C1, C2, X3N, Y3N)

  h1 = 2*pi/N; 
  h2 = (2*m*R+2*epsilon)/N; 
  alpha = k*cos(theta);
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up
  
  uij = inline('exp(1i*((alpha+j)*x1+sqrt(k^2-(alpha+j)^2)*x2)) + exp(1i*((alpha+j)*x1-sqrt(k^2-(alpha+j)^2)*x2))','x1','x2','k','alpha','j'); % incident field
  
  NearField = zeros(M1+M2+1);      
  index = Y(q_N~= 0)+N/2+N*(X(q_N~=0)+N/2 - 1);  

for j = -M1:1:M2         
    uij_ = uij(X*h1,Y*h2,k,alpha,j);  % grid values of incident field    
    uij_down = exp(1i*((alpha + j)*X(index)*h1 - sqrt(k^2 - (alpha + j)^2)*Y(index)*h2));   
    uij_up = exp(1i*((alpha + j)*X(index)*h1 + sqrt(k^2 - (alpha + j)^2)*Y(index)*h2));     
    [del1_usj del2_usj] = directSolver_GalerkinScheme_2D(N, k, m, R, theta, accuracy, q_3N, epsilon, B1, B2, C1, C2, X3N, Y3N, j);  
    
    for n = -M1:1:M2                                                                        
        int = q_N(index).*exp(-1i*((alpha + n).*X(index)*h1 + sqrt(k^2 - (alpha + n)^2).*(Y(index)*h2 - R))).*...
              ((alpha + n)./(4*pi*sqrt(k^2 - (alpha + n)^2)).*(del1_usj(index) + 1i*(alpha + j).*uij_(index))+...                 
              1/(4*pi).*(del2_usj(index) + 1i*sqrt(k^2 - (alpha + j)^2).*(uij_up - uij_down)));
                    
        NearField(n + M1 + 1,j + M1 + 1) = -h1*h2*sum(int(:));        
    end
end      

  w = zeros(1,M1+M2+1);         
  I = exp(-1i*R*sqrt(k^2 -(alpha + (-M1:M2)).^2));      
  w((alpha + (-M1:M2)).^2 < k^2) = 1i;  
  w((alpha + (-M1:M2)).^2 > k^2) = I((alpha + (-M1:M2)).^2 > k^2);     

  NearField = NearField.*(ones(M1+M2+1,1)*w.^-1).*(ones(M1+M2+1,1)*sqrt(k^2 -(alpha + (-M1:M2)).^2).^-1);
  
  