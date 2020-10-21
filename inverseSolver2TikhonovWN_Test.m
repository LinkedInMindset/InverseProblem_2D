

% This is the main program for the 2D periodic inverse solver including noise cases
% The parameters used for the papers are current in the module RUN_noise 

function  [picture q11_] = inverseSolver2TikhonovWN_Test(N, M1, M2, k, m , R, theta, accuracy, q11, q22, epsilon, structure, lastEigenvalues, noiselevel)


      h1 = 2*pi/N; 
      h2 = (2*m*R + 2*epsilon)/N;   
   alpha = k*cos(theta);    
    
   % Compute the coefficients for implementing VFFT and IVFFT
   B1=ones(N,1)*exp(-2*pi*1i*(1-N/2)*(1:N)/(N));
   B2=ones(N,1)*exp(-2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
   C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/(N));
   C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));    

   
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  
  % Compute the contrast q numerically
  q11_ = zeros(N); q22_ = zeros(N);
  if strcmp(structure,'sin')==1        
       q11_(abs(Y*h2) < -sin(X*h1)/2 + 1) = q11;   
       q22_(abs(Y*h2) < -sin(X*h1)/2 + 1) = q22;       
  elseif strcmp(structure,'cos')==1        
       q11_((Y*h2 < -cos(X*h1+1)/2 +1) & (Y*h2 > -cos(X*h1+1)/2 -0.5)) = q11; 
       q22_((Y*h2 < -cos(X*h1+1)/2 +1) & (Y*h2 > -cos(X*h1+1)/2 -0.5)) = q22;
  elseif strcmp(structure,'cube1')==1        
       q11_(((X*h1>=-pi & X*h1<=-pi/2) | (X*h1<=pi & X*h1>=pi/2)) & abs(Y*h2)<=1.5) = 0.5;
       q11_((X*h1>-pi/2 & X*h1<pi/2) & abs(Y*h2)<=0.5) = 1-3*1i;   
       q22_ = q11_;
  elseif strcmp(structure,'cube2')==1 
      %g = zeros(N);
      g = ((Y*h2 + 1).*(sin(X*h1).^2+0.5))/3 - 2*1i;
      
       q11_(((X*h1<-3*pi/4)|(X*h1>3*pi/4)) & abs(Y*h2) <=1.5) = g(((X*h1<-3*pi/4)|(X*h1>3*pi/4)) & abs(Y*h2) <=1.5);
       q11_((X*h1>=-pi/2 & X*h1<=pi/2) & abs(Y*h2) <=0.5) = g((X*h1>=-pi/2 & X*h1<=pi/2) & abs(Y*h2) <=0.5);        
       q11_((X*h1>=pi/2 & X*h1<=3*pi/4) & abs(Y*h2) <=4/pi*X*h1-1.5) = g((X*h1>=pi/2 & X*h1<=3*pi/4) & abs(Y*h2) <=4/pi*X*h1-1.5); 
       q11_((X*h1>=-3*pi/4 & X*h1<=-pi/2) & abs(Y*h2) <=-4/pi*X*h1-1.5) = g((X*h1>=-3*pi/4 & X*h1<=-pi/2) & abs(Y*h2) <=-4/pi*X*h1-1.5);       
       q22_ = q11_;
%        q_(((X*h1<-3*pi/4)|(X*h1>3*pi/4)) & abs(Y*h2) <=1.5) = q;
%        q_((X*h1>=-pi/2 & X*h1<=pi/2) & abs(Y*h2) <=0.5) = q;        
%        q_((X*h1>=pi/2 & X*h1<=3*pi/4) & abs(Y*h2) <=4/pi*X*h1-1.5) = q; 
%        q_((X*h1>=-3*pi/4 & X*h1<=-pi/2) & abs(Y*h2) <=-4/pi*X*h1-1.5) = q; 
  elseif strcmp(structure,'cube3')==1        
       q11_(((X*h1>=-pi & X*h1<=-pi/2) | (X*h1<=pi & X*h1>=pi/2)) & ((Y*h2<=1.5) & (Y*h2>=-0.5))) = 0.5;
       q11_((X*h1>-pi/2 & X*h1<pi/2) & ((Y*h2<=0.5) & (Y*h2>=-1.5))) = 1;   
       q22_ = q11_;
  elseif strcmp(structure,'cube4')==1        
       q11_(abs(X*h1)<=2 & Y*h2 <=1 & Y*h2 >=-1) = 0.5;       
       q22_ = q11_;  
  end    
  
  
    
   I = exp(-1i*R*sqrt(k^2 -(alpha + (-M1:M2)).^2));  
   w_star = zeros(1,M1+M2+1);   
   w_star((alpha + (-M1:M2)).^2 > k^2) = 1i;  
   w_star((alpha + (-M1:M2)).^2 < k^2) = I((alpha + (-M1:M2)).^2 < k^2);
      
  % Compute near field operator N from the Collocation direct solver  
  WNField = data_CollocationScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q11_, q22_, epsilon, B1, B2, C1, C2);
   
  WN = reshape(WNField(:,:,2,1),[M1+M2+1,M1+M2+1]);
  
  % Compute WN and WN_sharpe
  %   WN = -4*pi*NearField.*(ones(M1+M2+1,1)*w_star).';  
   ImWN = (WN - WN')./(2*1i);  
   ReWN = (WN + WN')/2;  
  [V,D] = eig(ReWN);  
  
  WN_shape = V*abs(D)*inv(V) - ImWN;   
  [V2, D2] = eig(WN_shape);  
  
  square_root_WN_shape = V2*sqrt(real(D2))*inv(V2);  
    
  
  % Compute artificial noise data  
   noise = rand(M1+M2+1);     
   square_root_WN_shapeNoise = square_root_WN_shape + noiselevel*noise*norm(square_root_WN_shape)/norm(noise) ; 
   [U,S,V] = svd(square_root_WN_shapeNoise);
   %norm(square_root_WN_shape - square_root_WN_shapeNoise)/norm(square_root_WN_shape)
          
   
  % Indice of the sampling domain
  index = Y(Y*h2<=R + 1 & Y*h2>=-R-1) + N/2 + N*(X(Y*h2<=R + 1 & Y*h2>=-R-1) + N/2 - 1);  
  
  % Exclude some lastEigenvalues may improve the picutre ? 
  index2 = (1 : 1 : M1+M2+1-lastEigenvalues);
 
  % Display the indice of propagating modes
  find(k^2 > (alpha + (-M1:M2)).^2) 
  
  % [sorted_D2 index3] = sort(sum(real(D2)),'descend'); % sort the eigenvalues in descending order
  % V2 = V2(:,index3); % eigenvectors in corresponding order    
  
  picture = ones(N)*10^40;
 
  % Compute regularization parameter gamma for some z
  r2 = 1i./(4*pi*sqrt(k^2 - (alpha + (-M1:M2)).^2)).*exp(-1i*((alpha + (-M1:M2)).*X(index(50))*h1 +...                 
                    sqrt(k^2 - (alpha + (-M1:M2)).^2).*(Y(index(50))*h2 - R)));                
   gamma = fzero(@(x) sum(abs(sum(conj(U).*(ones(M1+M2+1,1)*r2).')).^2./((sum(S).^2 + x).^2./(x^2 - (noiselevel)^2*sum(S).^2))), 0);
  
     for t = 1:length(index)                  
                                   
          r = 1i./(4*pi*sqrt(k^2 - (alpha + (-M1:M2)).^2)).*exp(-1i*((alpha + (-M1:M2)).*X(index(t))*h1 -...                 
                    sqrt(k^2 - (alpha + (-M1:M2)).^2).*(Y(index(t))*h2 + R)));   
                
      % Compute regularization parameter gamma for each z          
      % gamma = fzero(@(x)sum(abs(sum(conj(U).*(ones(M1+M2+1,1)*r).')).^2.*(x^2 - (noiselevel)^2*sum(S).^2)./((sum(S).^2 + x).^2)), 0); 
      
      picture(index(t)) = sum(abs(sum(conj(U(:,index2)).*(ones(length(index2),1)*r).')).^2./((sum(S(:,index2)).^2 + gamma).^2./sum(S(:,index2)).^2));     
             
     end
      

