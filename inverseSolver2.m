


function  [picture WN] = inverseSolver2(N, M1, M2, k, theta, structure, q, lastEigenvalues)


  R = 1.5; h1 = 2*pi/N; h2 = 4*R/N; alpha = k*cos(theta);   

   q_ = zeros(N);
    
  [X,Y] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  
  if structure == 1 
      
     q_(abs(Y*h2) < sin(2*X*h1)/2 + 1) = abs(X(abs(Y*h2) < sin(2*X*h1)/2 + 1)).^(1/2)*h1;
    
   elseif structure == 2
  
     q_(((X*h1>=-pi/4 & X*h1<pi/4) |(X*h1<-3*pi/4) | (X*h1>=3*pi/4)) & abs(Y*h2)<=1.5) = q;
  
     q_(((X*h1>=-3*pi/4 & X*h1<-pi/4) | (X*h1>=pi/4 & X*h1<3*pi/4)) & abs(Y*h2)<=0.5) = q;
   
   else

     q_(((X*h1>=-pi/8 & X*h1<=pi/8)|(X*h1<=-7*pi/8)|(X*h1>=7*pi/8)) & abs(Y*h2) <=1.5) = q;
     
     q_(((X*h1>=-3*pi/4 & X*h1<=-pi/4)|(X*h1>=pi/4 & X*h1<=3*pi/4)) & abs(Y*h2) <=0.5) = q;
 
     q_(((X*h1>=pi/8 & X*h1<=pi/4)) & abs(Y*h2) <=-8/pi*X*h1+2.5) = q;
  
     q_((X*h1>=-7*pi/8 & X*h1<=-3*pi/4) & abs(Y*h2) <=-8/pi*(X*h1+pi)+2.5) = q; 
 
     q_((X*h1>=-pi/4 & X*h1<=-pi/8) & abs(Y*h2) <=8/pi*X*h1+2.5) = q;

     q_((X*h1>=3*pi/4 & X*h1<=7*pi/8) & abs(Y*h2) <=8/pi*(X*h1-pi)+2.5) = q;
 
   end
        
  I = exp(-1i*R*sqrt(k^2 -(alpha + [-M1:M2]).^2));
  
  w_star = zeros(1,M1+M2+1); 
  
  w_star((alpha + [-M1:M2]).^2 > k^2) = 1i;
  
  w_star((alpha + [-M1:M2]).^2 < k^2) = I((alpha + [-M1:M2]).^2 < k^2);
  
  
  NearField = data2(N, M1, M2, k, theta, q_);       
  
  WN = -4*pi*NearField.*(ones(M1+M2+1,1)*w_star).';
  
  
  ImWN = (WN - WN')./(2*1i);
  
  ReWN = (WN + WN')/2;
  
  [V,D] = eig(ReWN);  
  
  WN_shape = V*abs(D)*inv(V) - ImWN; 
  
  [V2, D2] = eig(WN_shape); 
 
    
  index = Y(Y*h2<=R + 0.1 & Y*h2>=0) + N/2 + N*(X(Y*h2<=R + 0.1 & Y*h2>=0) + N/2 - 1);  % indice of tested domain
  
  index1 = find(k^2 > (alpha + [-M1:M2]).^2);  
 
  
  index1 - M1 - 1  % indice of propagating modes
  
  index2 = [1 : 1 : M1 + M2 + 1 - lastEigenvalues];
  
  [sorted_D2 index3] = sort(sum(real(D2)),'descend'); % sort the eigenvalues in descending order
  
  V2 = V2(:,index3); % eigenvectors in corresponding order
  
  
  picture = ones(N)*10^40;

     for t = 1:length(index)

          r = 1i./(4*pi*sqrt(k^2 - (alpha + [-M1:M2]).^2)).*exp(-1i*((alpha + [-M1:M2]).*X(index(t))*h1 +...
              sqrt(k^2 - (alpha + [-M1:M2]).^2).*(Y(index(t))*h2 - R)));  
          
           picture(index(t)) = sum(abs(sum(V2(:,index2).*(ones(length(index2),1)*r).')).^2./sorted_D2(:,index2));   
             
     end
    
   
    imagesc(-(-N/2+1:N/2)*h1,Y(Y(:,1)*h2<=2 & Y(:,1)*h2>=0)*h2,abs(picture(Y(:,1)*h2<=2 & Y(:,1)*h2>=0,1:N).^(-1/2)));    
    view(2); 
    set(gca,'YDir','normal')    
    axis tight
    colorbar
    shading interp    
     hold on
     
     t = linspace(-pi, 0,1000);
     ft = zeros(1,length(t));
     
     if structure ==1      
         
     ft = sin(2*t)./2 + 1;
     
     plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
     
     hold off
     
     elseif structure ==2
    
     ft((t>=-pi/4 & t<pi/4) |(t<-3*pi/4) | (t>=3*pi/4)) = 1.5;
     
     ft((t>=-3*pi/4 & t<-pi/4) | (t>=pi/4 & t<3*pi/4)) = 0.5;
     
     plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
     
     hold off
     
     else
         
     ft((t>=-pi/8 & t<=pi/8)|(t<=-7*pi/8)|(t>=7*pi/8)) = 1.5;

     ft((t>=-3*pi/4 & t<=-pi/4)|(t>=pi/4 & t<=3*pi/4)) = 0.5;
 
     ft(t>=pi/8 & t<=pi/4) =-8/pi*t(t>=pi/8 & t<=pi/4) + 2.5;
  
     ft(t>=-7*pi/8 & t<=-3*pi/4) = -8/pi*(t(t>=-7*pi/8 & t<=-3*pi/4)+pi) + 2.5; 
 
     ft(t>=-pi/4 & t<=-pi/8) = 8/pi*t(t>=-pi/4 & t<=-pi/8) + 2.5;

     ft(t>=3*pi/4 & t<=7*pi/8) = 8/pi*(t(t>=3*pi/4 & t<=7*pi/8)-pi) + 2.5;
     
     plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
     
     hold off
     
     end

    

    