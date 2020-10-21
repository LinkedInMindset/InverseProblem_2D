





function  picture = inverseSolver_CollocationScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q, epsilon, structure, lastEigenvalues)


  h1 = 2*pi/N; 
  h2 = (2*m*R + 2*epsilon)/N;   
  alpha = k*cos(theta);   
  [XN,YN] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
 
  B1=ones(N,1)*exp(-2*pi*1i*(1-N/2)*(1:N)/(N));
  B2=ones(N,1)*exp(-2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
  C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/(N));
  C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
   
  q_ = zeros(N);  
  if strcmp(structure,'sin')==1  
  q_(abs(YN*h2)<=sin(2*XN*h1)/2 + 0.5) = q;
  elseif strcmp(structure,'circle')==1   
  q_((XN*h1).^2+(YN*h2).^2 <= R^2) = q;
  elseif strcmp(structure,'rectangle1')==1 
  q_(abs(XN*h1)<=2 & abs(YN*h2)<=R) = q;   
  elseif strcmp(structure,'rectangle2')==1   
  %q_(((XN*h1>=-pi/4 & XN*h1<pi/4) |(XN*h1<-3*pi/4) | (XN*h1>=3*pi/4)) & abs(YN*h2)<=1.5) = q;  
  %q_(((XN*h1>=-3*pi/4 & XN*h1<-pi/4) | (XN*h1>=pi/4 & XN*h1<3*pi/4)) & abs(YN*h2)<=0.5) = q;    
  q_(((XN*h1>=-1.5 & XN*h1<1.5)) & abs(YN*h2)<=R) = q;  
  q_(((XN*h1<=-1.5) | (XN*h1>=1.5)) & abs(YN*h2)<=R/2) = q;    
  elseif strcmp(structure,'rectangle3')   
  q_(((XN*h1>=-pi/8 & XN*h1<=pi/8)|(XN*h1<=-7*pi/8)|(XN*h1>=7*pi/8)) & abs(YN*h2) <=1.5) = q;     
  q_(((XN*h1>=-3*pi/4 & XN*h1<=-pi/4)|(XN*h1>=pi/4 & XN*h1<=3*pi/4)) & abs(YN*h2) <=0.5) = q; 
  q_(((XN*h1>=pi/8 & XN*h1<=pi/4)) & abs(YN*h2) <=-8/pi*XN*h1+2.5) = q;  
  q_((XN*h1>=-7*pi/8 & XN*h1<=-3*pi/4) & abs(YN*h2) <=-8/pi*(XN*h1+pi)+2.5) = q;  
  q_((XN*h1>=-pi/4 & XN*h1<=-pi/8) & abs(YN*h2) <=8/pi*XN*h1+2.5) = q;
  q_((XN*h1>=3*pi/4 & XN*h1<=7*pi/8) & abs(YN*h2) <=8/pi*(XN*h1-pi)+2.5) = q; 
  end
  
  I = exp(-1i*R*sqrt(k^2-(alpha+(-M1:M2)).^2));  
  w_star = zeros(1,M1+M2+1);   
  w_star((alpha + (-M1:M2)).^2 > k^2) = 1i;  
  w_star((alpha + (-M1:M2)).^2 < k^2) = I((alpha + (-M1:M2)).^2 < k^2);
  
  
  NearField = data_CollocationScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q_, epsilon, B1, B2, C1, C2);   
  WN = -4*pi*NearField.*(ones(M1+M2+1,1)*w_star).';
  
  
  ImWN = (WN - WN')./(2*1i);  
  ReWN = (WN + WN')/2;  
  [V,D] = eig(ReWN);    
  WN_shape = V*abs(D)*inv(V) - ImWN;   
  [V2, D2] = eig(WN_shape); 
 
    
  index = YN(YN*h2<=R + 0.5 & YN*h2>=0) + N/2 + N*(XN(YN*h2<=R + 0.5 & YN*h2>=0) + N/2 - 1);  % indice of tested domain
  
  index1 = find(k^2 > (alpha + (-M1:M2)).^2);    
  index1 - M1 - 1  % indice of propagating modes
  
  index2 = (1 : 1 : M1 + M2 + 1 - lastEigenvalues);
  
  [sorted_D2 index3] = sort(sum(real(D2)),'descend'); % sort the eigenvalues in descending order  
  V2 = V2(:,index3); % eigenvectors in corresponding order  
  
  picture = ones(N)*10^40;
     for t = 1:length(index)
          r = 1i./(4*pi*sqrt(k^2 - (alpha + (-M1:M2)).^2)).*exp(-1i*((alpha + (-M1:M2)).*XN(index(t))*h1 +...
              sqrt(k^2 - (alpha + (-M1:M2)).^2).*(YN(index(t))*h2 - R)));  
          
           picture(index(t)) = sum(abs(sum(V2(:,index2).*(ones(length(index2),1)*r).')).^2./sorted_D2(:,index2));               
     end
    
   
    imagesc(-(-N/2+1:N/2)*h1,YN(YN(:,1)*h2<=2*pi & YN(:,1)*h2>=0)*h2,abs(picture(YN(:,1)*h2<=2*pi & YN(:,1)*h2>=0,1:N).^(-1/2)));    
    view(2); 
    set(gca,'YDir','normal')    
    axis tight
    colorbar
    shading interp    
    hold on
     
%      t = linspace(-pi, 0,1000);
%      ft = zeros(1,length(t));       
%      ft = sin(2*t)./2 + 0.5;      
%      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
     
    