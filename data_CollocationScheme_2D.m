



function   [WNearField,w_starPlus,w_starMinus] = data_CollocationScheme_2D(N, M1, M2, k, m , R, R_meas, theta, accuracy, q11_, q22_, epsilon, B1, B2, C1, C2)

     h1 = 2*pi/N; 
     h2 = (2*m*R+2*epsilon)/N;   
  alpha = k*cos(theta);
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  

  
  WNearField = zeros(M1+M2+1,M1+M2+1,2,2);    
  index = Y(q11_~= 0)+N/2+N*(X(q11_~=0)+N/2-1);  
  
  
  I = exp(-1i*R_meas*sqrt(k^2 -(alpha + (-M1:M2)).^2));      
  
  w_starPlus = zeros(1,M1+M2+1);        
  w_starPlus((alpha + (-M1:M2)).^2 > k^2) = 1i;  
  w_starPlus((alpha + (-M1:M2)).^2 < k^2) = I((alpha + (-M1:M2)).^2 < k^2); 
  
  w_starMinus = zeros(1,M1+M2+1);    
  w_starMinus((alpha + (-M1:M2)).^2 > k^2) = 1i;  
  w_starMinus((alpha + (-M1:M2)).^2 < k^2) = 1i*I((alpha + (-M1:M2)).^2 < k^2); 
  
    
  w_Plus = zeros(1,M1+M2+1);        
  w_Plus((alpha + (-M1:M2)).^2 < k^2) = 1i;  
  w_Plus((alpha + (-M1:M2)).^2 > k^2) = I((alpha + (-M1:M2)).^2 > k^2); 
  
  w_Minus = zeros(1,M1+M2+1);    
  w_Minus((alpha + (-M1:M2)).^2 < k^2) = 1;  
  w_Minus((alpha + (-M1:M2)).^2 > k^2) = I((alpha + (-M1:M2)).^2 > k^2); 
  
   A_Plus = w_Plus.^-1.*sqrt(k^2-(alpha + (-M1:M2)).^2).^-1;
  A_Minus = w_Minus.^-1.*sqrt(k^2-(alpha + (-M1:M2)).^2).^-1;
  
  
for incField = 1:2
    
    if incField == 1
        A = A_Plus;
    else
        A = A_Minus;
    end

    for j = -M1:1:M2         
    
    [del1_utj,del2_utj] = directSolver_CollocationScheme_2D(N, k, m, R, theta, accuracy, q11_, q22_, epsilon, B1, B2, C1, C2, incField, j);     
      
    for n = -M1:1:M2       
                        
            alphan = alpha + n;
             betan = sqrt(k^2 - alphan^2); 
                
           rn_plus = 1i/(4*pi*betan)*exp(-1i*(alphan*X*h1 + betan*(Y*h2-R_meas)));        
          rn_minus = 1i/(4*pi*betan)*exp(-1i*(alphan*X*h1 - betan*(Y*h2+R_meas)));  
        
          int_plus = -1i*alphan.*rn_plus(index).*q11_(index).*del1_utj(index) + -1i*betan.*rn_plus(index).*q22_(index).*del2_utj(index);
         int_minus = -1i*alphan.*rn_minus(index).*q11_(index).*del1_utj(index) + 1i*betan.*rn_minus(index).*q22_(index).*del2_utj(index);
                    
        %WNearField(n + M1 + 1,j + M1 + 1,1,incField) = -h1*h2*w_starPlus(n + M1 + 1).*A(j + M1 + 1).*(sum(int_plus(:))); 
        %WNearField(n + M1 + 1,j + M1 + 1,2,incField) = -h1*h2*w_starPlus(n + M1 + 1).*A(j + M1 + 1).*(sum(int_plus(:))); 
        WNearField(n + M1 + 1,j + M1 + 1,1,incField) = h1*h2*w_starPlus(n + M1 + 1).*A(j + M1 + 1).*(sum(int_plus(:)) + sum(int_minus(:)));   
        WNearField(n + M1 + 1,j + M1 + 1,2,incField) = h1*h2*w_starMinus(n + M1 + 1).*A(j + M1 + 1).*(sum(int_plus(:)) - sum(int_minus(:)));    
    end
    end      
end

  WNearField = 4*pi*WNearField; 
  %.*(ones(M1+M2+1,1)*w.^-1).*(ones(M1+M2+1,1)*sqrt(k^2 -(alpha + (-M1:M2)).^2).^-1);
  
  
  
  
  
  
  