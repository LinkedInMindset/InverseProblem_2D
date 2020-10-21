
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to image periodic structures using Factorization method. 
% N is parameter of discretization
% M1 and M2 are lower and upper bounds of for the number of the incident fields ui_j that means j varies from M1 to M2. 
% M1 + M2 is the number of the incident fields used
% k is wave number
% [-pi,pi]*[-m*R-epsilon,m*R+epsilon] is one period. A good example: m = 2.5, R = 1, epsilon = 0.5
% structure must be bounded by [-R,R] in x2-direction 
% theta is angle of incident plane wave ui_j
% alpha = k*cos(theta)
% plane wave ui_j = exp(1i*((alpha+j)*x1+sqrt(k^2-(alpha+j)^2)*x2)) + exp(1i*((alpha+j)*x1-sqrt(k^2-(alpha+j)^2)*x2))
% accuracy is the residual of GMRES iteration
% q is the contrast 
% [XN,YN]=meshgrid(-N/2+1:N/2)
% q_3N and q_N are the grid values the contrast q on grid 3N and N, respectively 
% structure is the parameter for kind of periodic structure. It must be an element of {'sin','circle','rectangle1','rectangle2'} 
% lastEigenvalues is the number of last eigenvalues that we want to exclude. That's just a trick to try to improve the quality of the rescontructed image
% lastEigenvalues is 0 that means you retains all eigenvalues  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function  picture = inverseSolver_GalerkinScheme_2D(N, M1, M2, k, m , R, theta, accuracy, q, epsilon, structure, lastEigenvalues)


   h1 = 2*pi/N; 
   h2 = (2*m*R + 2*epsilon)/N;   
   alpha = k*cos(theta);   
   
   B1=ones(3*N,1)*exp(-2*pi*1i*(1-3*N/2)*(1:3*N)/(3*N));
   B2=ones(3*N,1)*exp(-2*pi*1i*(((1-3*N/2)*(1:3*N)+(3*N)^2/4-1))/(3*N));
   C1=ones(3*N,1)*exp(2*pi*1i*(1-3*N/2)*(1:3*N)/(3*N));
   C2=ones(3*N,1)*exp(2*pi*1i*(((1-3*N/2)*(1:3*N)+(3*N)^2/4-1))/(3*N));     
    
  [XN,YN] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  [X3N,Y3N] = meshgrid(-3*N/2+1:3*N/2);  
  
  q_N = zeros(N);
  q_3N = zeros(3*N);    
  
  if strcmp(structure,'sin')==1  
  q_3N(abs(Y3N*h2/3)<=sin(2*X3N*h1/3)/2 + 0.5) = q;
  q_N(abs(YN*h2)<=sin(2*XN*h1)/2 + 0.5) = q;
  elseif strcmp(structure,'circle')==1   
  q_3N((X3N*h1/3).^2+(Y3N*h2/3).^2 <= R^2) = q;
  q_N((XN*h1).^2+(YN*h2).^2 <= R^2) = q;
  elseif strcmp(structure,'rectangle1')==1 
  q_3N(abs(X3N*h1/3)<=2 & abs(Y3N*h2/3)<=R) = q;   
  q_N(abs(XN*h1)<=2 & abs(YN*h2)<=R) = q; 
  elseif strcmp(structure,'rectangle2')==1   
  %q_(((XN*h1>=-pi/4 & XN*h1<pi/4) |(XN*h1<-3*pi/4) | (XN*h1>=3*pi/4)) & abs(YN*h2)<=1.5) = q;  
  %q_(((XN*h1>=-3*pi/4 & XN*h1<-pi/4) | (XN*h1>=pi/4 & XN*h1<3*pi/4)) & abs(YN*h2)<=0.5) = q;    
  q_3N(((X3N*h1/3>=-1.5 & X3N*h1/3<1.5)) & abs(Y3N*h2/3)<=R) = q;  
  q_3N(((X3N*h1/3<=-1.5) | (X3N*h1/3>=1.5)) & abs(Y3N*h2/3)<=R/2) = q;    
  q_N(((XN*h1>=-1.5 & XN*h1<1.5)) & abs(YN*h2)<=R) = q;  
  q_N(((XN*h1<=-1.5) | (XN*h1>=1.5)) & abs(YN*h2)<=R/2) = q;   
%   elseif strcmp(structure,'rectangle3')   
%   q_3N(((X3N*h1/3>=-pi/8 & X3N*h1/3<=pi/8)|(X3N*h1/3<=-7*pi/8)|(X3N*h1/3>=7*pi/8)) & abs(Y3N*h2/3) <=1.5) = q;     
%   q_3N(((X3N*h1/3>=-3*pi/4 & X3N*h1/3<=-pi/4)|(X3N*h1/3>=pi/4 & X3N*h1/3<=3*pi/4)) & abs(Y3N*h2/3) <=0.5) = q; 
%   q_3N(((X3N*h1/3>=pi/8 & X3N*h1/3<=pi/4)) & abs(Y3N*h2/3) <=-8/pi*X3N*h1/3+2.5) = q;  
%   q_3N((X3N*h1/3>=-7*pi/8 & X3N*h1/3<=-3*pi/4) & abs(Y3N*h2/3) <=-8/pi*(X3N*h1/3+pi)+2.5) = q;  
%   q_3N((X3N*h1/3>=-pi/4 & X3N*h1/3<=-pi/8) & abs(Y3N*h2/3) <=8/pi*X3N*h1/3+2.5) = q;
%   q_3N((X3N*h1/3>=3*pi/4 & X3N*h1/3<=7*pi/8) & abs(Y3N*h2/3) <=8/pi*(X3N*h1/3-pi)+2.5) = q; 
%   
%   q_N(((XN*h1>=-pi/8 & XN*h1<=pi/8)|(XN*h1<=-7*pi/8)|(XN*h1>=7*pi/8)) & abs(YN*h2) <=1.5) = q;     
%   q_N(((XN*h1>=-3*pi/4 & XN*h1<=-pi/4)|(XN*h1>=pi/4 & XN*h1<=3*pi/4)) & abs(YN*h2) <=0.5) = q; 
%   q_N(((XN*h1>=pi/8 & XN*h1<=pi/4)) & abs(YN*h2) <=-8/pi*XN*h1+2.5) = q;  
%   q_N((XN*h1>=-7*pi/8 & XN*h1<=-3*pi/4) & abs(YN*h2) <=-8/pi*(XN*h1+pi)+2.5) = q;  
%   q_N((XN*h1>=-pi/4 & XN*h1<=-pi/8) & abs(YN*h2) <=8/pi*XN*h1+2.5) = q;
%   q_N((XN*h1>=3*pi/4 & XN*h1<=7*pi/8) & abs(YN*h2) <=8/pi*(XN*h1-pi)+2.5) = q; 
  end
  
        
  I = exp(-1i*R*sqrt(k^2 -(alpha + (-M1:M2)).^2));  
  w_star = zeros(1,M1+M2+1);   
  w_star((alpha + (-M1:M2)).^2 > k^2) = 1i;  
  w_star((alpha + (-M1:M2)).^2 < k^2) = I((alpha + (-M1:M2)).^2 < k^2);
  
  
  NearField = data_GalerkinScheme_2D(N, M1, M2, k, m, R, theta, accuracy, q_N, q_3N, epsilon, B1, B2, C1, C2, X3N, Y3N);   
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
    %view(2); 
    set(gca,'YDir','normal')    
    axis tight
    colorbar
    shading interp    
     
     
%      t = linspace(-pi, 0,1000);
%      ft = zeros(1,length(t));
%      
%      if structure ==1      
%          
%      ft = sin(2*t)./2 + 1;
%      
%      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
%      
%      hold off
%      
%      elseif structure ==2
%     
%      ft((t>=-pi/4 & t<pi/4) |(t<-3*pi/4) | (t>=3*pi/4)) = 1.5;
%      
%      ft((t>=-3*pi/4 & t<-pi/4) | (t>=pi/4 & t<3*pi/4)) = 0.5;
%      
%      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
%      
%      hold off
%      
%      else
%          
%      ft((t>=-pi/8 & t<=pi/8)|(t<=-7*pi/8)|(t>=7*pi/8)) = 1.5;
% 
%      ft((t>=-3*pi/4 & t<=-pi/4)|(t>=pi/4 & t<=3*pi/4)) = 0.5;
%  
%      ft(t>=pi/8 & t<=pi/4) =-8/pi*t(t>=pi/8 & t<=pi/4) + 2.5;
%   
%      ft(t>=-7*pi/8 & t<=-3*pi/4) = -8/pi*(t(t>=-7*pi/8 & t<=-3*pi/4)+pi) + 2.5; 
%  
%      ft(t>=-pi/4 & t<=-pi/8) = 8/pi*t(t>=-pi/4 & t<=-pi/8) + 2.5;
% 
%      ft(t>=3*pi/4 & t<=7*pi/8) = 8/pi*(t(t>=3*pi/4 & t<=7*pi/8)-pi) + 2.5;
%      
%      plot(t,ft,'--','LineWidth',2,'Color',[.9 0.9 0.9]);    
%      
%      hold off
%      
%      end