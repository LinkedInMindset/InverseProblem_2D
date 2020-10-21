

function  picture = inverseSolver2TikhonovWN(N, M1, M2, k, m , R, theta, accuracy, q, epsilon, structure, propagateModes, noiselevel)


      h1 = 2*pi/N; 
      h2 = (2*m*R + 2*epsilon)/N;   
   alpha = k*cos(theta);    
    
   B1=ones(N,1)*exp(-2*pi*1i*(1-N/2)*(1:N)/(N));
   B2=ones(N,1)*exp(-2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
   C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/(N));
   C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));    

  q_ = zeros(N);
  
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  
  if strcmp(structure,'sin')==1        
       q_(abs(Y*h2) < -sin(X*h1)/2 +1) = q;
    
  elseif strcmp(structure,'cube1')==1        
       q_(((X*h1>=-pi & X*h1<=-pi/2) | (X*h1<=pi & X*h1>=pi/2)) & abs(Y*h2)<=1.5) = q;
       q_((X*h1>-pi/2 & X*h1<pi/2) & abs(Y*h2)<=0.5) = q;
  
%      q_(((X*h1>=-pi/4 & X*h1<pi/4) |(X*h1<-3*pi/4) | (X*h1>=3*pi/4)) & abs(Y*h2)<=1.5) = q;  
%      q_(((X*h1>=-3*pi/4 & X*h1<-pi/4) | (X*h1>=pi/4 & X*h1<3*pi/4)) & abs(Y*h2)<=0.5) = q;
   
   elseif strcmp(structure,'cube2')==1 
        q_(((X*h1<-3*pi/4)|(X*h1>3*pi/4)) & abs(Y*h2) <=1.5) = q;
        q_((X*h1>=-pi/2 & X*h1<=pi/2) & abs(Y*h2) <=0.5) = q;        
        q_((X*h1>=pi/2 & X*h1<=3*pi/4) & abs(Y*h2) <=4/pi*X*h1-1.5) = q; 
        q_((X*h1>=-3*pi/4 & X*h1<=-pi/2) & abs(Y*h2) <=-4/pi*X*h1-1.5) = q; 
        
%        q_(((X*h1>=-pi/8 & X*h1<=pi/8)|(X*h1<-7*pi/8)|(X*h1>7*pi/8)) & abs(Y*h2) <=1.5) = q;
%        q_(((X*h1>=-3*pi/4 & X*h1<-pi/4)|(X*h1>=pi/4 & X*h1<3*pi/4)) & abs(Y*h2) <=0.5) = q; 
%        q_(((X*h1>pi/8 & X*h1<=pi/4)) & abs(Y*h2) <=-8/pi*X*h1+2.5) = q;  
%        q_((X*h1>=-7*pi/8 & X*h1<-3*pi/4) & abs(Y*h2) <=-8/pi*(X*h1+pi)+2.5) = q;  
%        q_((X*h1>=-pi/4 & X*h1<-pi/8) & abs(Y*h2) <=8/pi*X*h1+2.5) = q;
%        q_((X*h1>=3*pi/4 & X*h1<=7*pi/8) & abs(Y*h2) <=8/pi*(X*h1-pi)+2.5) = q; 
  end    
    
   I = exp(-1i*R*sqrt(k^2 -(alpha + (-M1:M2)).^2));  
   w_star = zeros(1,M1+M2+1);   
   w_star((alpha + (-M1:M2)).^2 > k^2) = 1i;  
   w_star((alpha + (-M1:M2)).^2 < k^2) = I((alpha + (-M1:M2)).^2 < k^2);
      
  % Data generated from direct problems
  
  NearField = data_CollocationScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q_, epsilon, B1, B2, C1, C2);   
  
     WN = -4*pi*NearField.*(ones(M1+M2+1,1)*w_star).';  
   ImWN = (WN - WN')./(2*1i);  
   ReWN = (WN + WN')/2;  
  [V,D] = eig(ReWN);  
  
  WN_shape = V*abs(D)*inv(V) - ImWN;   
  [V2, D2] = eig(WN_shape);  
  
  square_root_WN_shape = V2*sqrt(real(D2))*inv(V2);  
    
  
  % Create artificial noise data
  
   noise = rand(M1+M2+1);  
   
   square_root_WN_shapeNoise = square_root_WN_shape + noiselevel*noise*norm(square_root_WN_shape)/norm(noise) ;
 
  [VNoise3 DNoise3] = eig(square_root_WN_shapeNoise*square_root_WN_shapeNoise');
        
  index = Y(Y*h2<=R + 0.1 & Y*h2>=0) + N/2 + N*(X(Y*h2<=R + 0.1 & Y*h2>=0) + N/2 - 1);  % indice of tested domain
  
  if strcmp(propagateModes,'yes') ==1  
  index2 = find(k^2 > (alpha + (-M1:M2)).^2);  % indice of propagating modes  
  elseif strcmp(propagateModes,'no') ==1
  index2 = (1 : 1 : M1+M2+1);
  end
    
%   [sorted_D2 index3] = sort(sum(real(D2)),'descend'); % sort the eigenvalues in descending order
%   V2 = V2(:,index3); % eigenvectors in corresponding order    
  
  picture = ones(N)*10^40;

  r2 = 1i./(4*pi*sqrt(k^2 - (alpha + (-M1:M2)).^2)).*exp(-1i*((alpha + (-M1:M2)).*X(index(5))*h1 +...                 
                    sqrt(k^2 - (alpha + (-M1:M2)).^2).*(Y(index(5))*h2 - R)));
                
   gamma = fzero(@(x) sum(abs(sum(VNoise3.*conj((ones(M1+M2+1,1)*r2)).')).^2./((sum(DNoise3) + x).^2./(x^2 - (noiselevel)^2*sum(DNoise3)))), 0);
  
     for t = 1:length(index)         

          r = 1i./(4*pi*sqrt(k^2 - (alpha + (-M1:M2)).^2)).*exp(-1i*((alpha + (-M1:M2)).*X(index(t))*h1 +...                 
                    sqrt(k^2 - (alpha + (-M1:M2)).^2).*(Y(index(t))*h2 - R)));                                       

          picture(index(t)) = sum(abs(sum(VNoise3(:,index2).*conj((ones(length(index2),1)*r)).')).^2./((sum(DNoise3(:,index2)) + gamma).^2./sum(DNoise3(:,index2))));     
             
     end
      
   sorted_DNoise3 = sort(sum(DNoise3),'descend');   
      picture_Per = cat(2,picture,picture);
     
    %subplot(121)
    imagesc((-N/2+1:3*N/2)*h1,Y(Y(:,1)*h2<=2 & Y(:,1)*h2>=0)*h2,abs(picture_Per(Y(:,1)*h2<=2 & Y(:,1)*h2>=0,:).^(-1/3)));    
    view(2); 
    set(gca,'YDir','normal')   
    axis tight
    colorbar
    shading interp    
    hold on
    
    t = linspace(-pi, pi,1000);
     ft = zeros(1,length(t));
     
     if strcmp(structure,'sin')==1               
      ft = -sin(t)/2 + 1;     
      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);         
      hold off
     
     elseif strcmp(structure,'cube1')==1          
      ft((t>=-pi & t<=-pi/2) | (t<=pi & t>=pi/2)) = 1.5;
      ft(t>-pi/2 & t<pi/2) = 0.5;     
%     ft((t>=-pi/4 & t<pi/4) |(t<-3*pi/4) | (t>=3*pi/4)) = 1.5;     
%     ft((t>=-3*pi/4 & t<-pi/4) | (t>=pi/4 & t<3*pi/4)) = 0.5;     
      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);        
      hold off
     
     elseif strcmp(structure,'cube2')==1      
                  
        ft((t<-3*pi/4)|(t>3*pi/4)) = 1.5;
        ft(t>=-pi/2 & t<=pi/2) = 0.5;        
        ft(t>=pi/2 & t<=3*pi/4)= 4/pi*t(t>=pi/2 & t<=3*pi/4)-1.5; 
        ft(t>=-3*pi/4 & t<=-pi/2)= -4/pi*t(t>=-3*pi/4 & t<=-pi/2)-1.5; 
        
%       ft((t>=-pi/8 & t<=pi/8)|(t<-7*pi/8)|(t>7*pi/8)) = 1.5;
%       ft((t>=-3*pi/4 & t<-pi/4)|(t>=pi/4 & t<3*pi/4)) = 0.5; 
%       ft(t>pi/8 & t<=pi/4) = -8/pi*t(t>pi/8 & t<=pi/4) + 2.5;  
%       ft(t>=-7*pi/8 & t<-3*pi/4) = -8/pi*(t(t>=-7*pi/8 & t<-3*pi/4)+pi) + 2.5;  
%       ft(t>=-pi/4 & t<-pi/8) = 8/pi*t(t>=-pi/4 & t<-pi/8) + 2.5;
%       ft(t>=3*pi/4 & t<=7*pi/8) = 8/pi*(t(t>=3*pi/4 & t<=7*pi/8)-pi) + 2.5;     
      plot(t, ft,'--','LineWidth',2 ,'color', [.9 0.9 0.9]);         
      hold off
     
     end
   % hold on    
   % semilogy([-M1:M2],abs(sorted_DNoise3),'bo')    
   % hold off