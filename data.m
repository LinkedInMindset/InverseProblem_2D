


function   NearField = data(N, M1, M2, k, theta, structure, q)



  R = 1.5;  h1 = 2*pi/N; h2 = 4*R/N;   

  ui = inline('exp(1i*((k*cos(theta)+j)*x1+sqrt(k^2-(k*cos(theta)+j)^2)*x2)) + exp(1i*((k*cos(theta)+j)*x1-sqrt(k^2-(k*cos(theta)+j)^2)*x2))','x1','x2','k','theta','j'); % incident field

  u = zeros(N); 
  
  NearField = zeros(M1+M2+1);  

  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up
  
    q_ = zeros(N);
    
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  
  if structure == 1 
      
     q_(abs(Y*h2) < sin(2*X*h1)/2 + 1) = q;
    
   elseif structure == 2
  
     q_(((X*h1>=-pi/4 & X*h1<pi/4) |(X*h1<-3*pi/4) | (X*h1>=3*pi/4)) & abs(Y*h2)<=1.5) = q;
  
     q_(((X*h1>=-3*pi/4 & X*h1<-pi/4) | (X*h1>=pi/4 & X*h1<3*pi/4)) & abs(Y*h2)<=0.5) = q;
   
   else

     q_(((X*h1>=-pi/8 & X*h1<=pi/8)|(X*h1<-7*pi/8)|(X*h1>7*pi/8)) & abs(Y*h2) <=1.5) = q;

     q_(((X*h1>=-3*pi/4 & X*h1<-pi/4)|(X*h1>=pi/4 & X*h1<3*pi/4)) & abs(Y*h2) <=0.5) = q;
 
     q_(((X*h1>pi/8 & X*h1<=pi/4)) & abs(Y*h2) <=-8/pi*X*h1+2.5) = q;
  
     q_((X*h1>=-7*pi/8 & X*h1<-3*pi/4) & abs(Y*h2) <=-8/pi*(X*h1+pi)+2.5) = q; 
 
     q_((X*h1>=-pi/4 & X*h1<-pi/8) & abs(Y*h2) <=8/pi*X*h1+2.5) = q;

     q_((X*h1>=3*pi/4 & X*h1<=7*pi/8) & abs(Y*h2) <=8/pi*(X*h1-pi)+2.5) = q;
 
   end
  
  alpha = k*cos(theta);

  %index = Y(Y*h2<R & Y*h2>=-R) + N/2 + N*(X(Y*h2<R & Y*h2>=-R) + N/2 - 1);  
  
  index = Y(q_ ~= 0) + N/2 + N*(X(q_ ~=0) + N/2 - 1);  

for j = -M1:1:M2
    
    [del1u del2u] = directSolverH2(N, k, theta, j, q_);
    
   
    down = exp(1i*((alpha + j)*X(index)*h1 - sqrt(k^2 - (alpha + j)^2)*Y(index)*h2));
    
  
    up = exp(1i*((alpha + j)*X(index)*h1 + sqrt(k^2 - (alpha + j)^2)*Y(index)*h2));
    
   
    ui_ = ui(X*h1, Y*h2, k, theta, j);  % grid values of incident field
  
    for n = -M1:1:M2     
                                                                   
        int = q_(index).*exp(-1i*((alpha + n).*X(index)*h1 + sqrt(k^2 - (alpha + n)^2).*(Y(index)*h2 - R))).*...
              ((alpha + n)./(4*pi*sqrt(k^2 - (alpha + n)^2)).*(del1u(index) + 1i*(alpha + j).*ui_(index))+...                 
              1/(4*pi).*(del2u(index) + 1i*sqrt(k^2 - (alpha + j)^2).*(up - down)));
                    
        NearField(n + M1 + 1,j + M1 + 1) = -h1*h2*sum(int(:));
        
    end
end      

  w = zeros(1,M1+M2+1);    
     
  I = exp(-1i*R*sqrt(k^2 -(alpha + [-M1:M2]).^2));  
    
  w((alpha + [-M1:M2]).^2 < k^2) = 1i;
  
  w((alpha + [-M1:M2]).^2 > k^2) = I((alpha + [-M1:M2]).^2 > k^2);     

  NearField = NearField.*(ones(M1+M2+1,1)*w.^-1).*(ones(M1+M2+1,1)*sqrt(k^2 -(alpha + [-M1:M2]).^2).^-1);
  
  