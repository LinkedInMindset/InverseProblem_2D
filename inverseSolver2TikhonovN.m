

function  picture = inverseSolver2TikhonovN(N, M1, M2, k, theta, structure, q, lastEigenvalues, noiselevel)


  R = 1.5;  h1 = 2*pi/N; h2 = 4*R/N; alpha = k*cos(theta);
  

  w_star = zeros(1,M1+M2+1);  q_ = zeros(N);
  
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
    
  I = exp(-1i*R*sqrt(k^2 -(alpha + [-M1:M2]).^2));
  
  
  w_star((alpha + [-M1:M2]).^2 > k^2) = 1i;
  
  w_star((alpha + [-M1:M2]).^2 < k^2) = I((alpha + [-M1:M2]).^2 < k^2);
  
  
  % Data  
  
  NearField = data2(N, M1, M2, k, theta, q_);   
  
  WN = -4*pi*NearField.*(ones(M1+M2+1,1)*w_star).';
  
  ImWN = (WN - WN')./(2*1i);
  
  ReWN = (WN + WN')/2;
  
  [V,D] = eig(ReWN);  
  
  WN_shape = V*abs(D)*inv(V) - ImWN; 
  
  [V2, D2] = eig(WN_shape);  
  
  square_root_WN_shape = V2*sqrt(real(D2))*inv(V2);
  
    
  
  % Create artificial noise data
  
   noise = rand(M1+M2+1);
  
  NearFieldNoise = NearField + noiselevel*noise*norm(NearField)/norm(noise);       
    
  WNNoise = -4*pi*NearFieldNoise.*(ones(M1+M2+1,1)*w_star).';
  
  
  ImWNNoise = (WNNoise - WNNoise')./(2*1i);
  
  ReWNNoise = (WNNoise + WNNoise')/2;
  
  [VNoise,DNoise] = eig(ReWNNoise);  
  
  WN_shapeNoise = VNoise*abs(DNoise)*inv(VNoise) - ImWNNoise; 
  
  [VNoise2, DNoise2] = eig(WN_shapeNoise);
  
  [sorted_D2 index3] = sort(sum(abs(DNoise2)),'descend'); % sort the eigenvalues in descending order
  
  square_root_WN_shapeNoise = VNoise2*sqrt(real(DNoise2))*inv(VNoise2);
  
  
  noiselevel2 = norm(square_root_WN_shapeNoise - square_root_WN_shape)/norm(square_root_WN_shape)
  
  
  [VNoise3 DNoise3] = eig(square_root_WN_shapeNoise*square_root_WN_shapeNoise');
  
    
  index = Y(Y*h2<=R + 0.5 & Y*h2>=0) + N/2 + N*(X(Y*h2<=R + 0.5 & Y*h2>=0) + N/2 - 1);  % indice of tested domain
  
 % index1 = find(k^2 > (alpha + [-M1:M2]).^2);  % indice of propagating modes
  
  index2 = [1 : 1 : M1+M2+1 - lastEigenvalues];
  
%   [sorted_D2 index3] = sort(sum(real(D2)),'descend'); % sort the eigenvalues in descending order
%   
%   V2 = V2(:,index3); % eigenvectors in corresponding order  
  
  
  picture = ones(N)*10^40;

  r2 = 1i./(4*pi*sqrt(k^2 - (alpha + [-M1:M2]).^2)).*exp(-1i*((alpha + [-M1:M2]).*X(index(8))*h1 +...                 
                    sqrt(k^2 - (alpha + [-M1:M2]).^2).*(Y(index(8))*h2 - R)));
                
   gamma = fzero(@(x) sum(abs(sum(VNoise3.*(ones(M1+M2+1,1)*r2).')).^2./((sum(DNoise3) + x).^2./(x^2 - (noiselevel2)^2*sum(DNoise3)))), 0);
  
     for t = 1:length(index)
         
%           r = 1i./(4*pi*sqrt(k^2 - (alpha + [-M/2 + 1:M/2]).^2)).*exp(-1i*((alpha + [-M/2 + 1:M/2]).*X(index(t))*h1 +...
%               sqrt(k^2 - (alpha + [-M/2 + 1:M/2]).^2).*(Y(index(t))*h2 - R)));
%           
             
          
%    Z(index(t)) = sum(abs(sum(V3.*(ones(M,1)*r).')).^2./((sum(D3) + gamma).^2./sum(D3)));    
         

          r = 1i./(4*pi*sqrt(k^2 - (alpha + [-M1:M2]).^2)).*exp(-1i*((alpha + [-M1:M2]).*X(index(t))*h1 +...                 
                    sqrt(k^2 - (alpha + [-M1:M2]).^2).*(Y(index(t))*h2 - R)));
                
          
                

          picture(index(t)) = sum(abs(sum(VNoise3(:,index2).*(ones(length(index2),1)*r).')).^2./((sum(DNoise3(:,index2)) + gamma).^2./sum(DNoise3(:,index2))));     
             
     end
     
  
    
   
  % sorted_DNoise3 = sort(sum(DNoise3),'descend');
     
    subplot(121)
    imagesc(-(-N/2+1:N/2)*h1,Y(Y(:,1)*h2<=3 & Y(:,1)*h2>=0)*h2,abs(picture(Y(:,1)*h2<=3 & Y(:,1)*h2>=0,1:N).^(-1/2)));    
    view(2); 
    set(gca,'YDir','normal')
    axis square
    axis tight
    colorbar
    shading interp    
    hold on
    
    t = linspace(-pi, pi,1000); 
    ft = sin(2*t)./2 + 1;
    plot(t,ft, 'LineWidth',2,'Color',[.6 0 0]);    
    hold off
    
    subplot(122)     
    
    semilogy([-M1:M2],abs(sorted_D2),'ro')
  %  hold on
   % subplot(122)
    %hold on
    
   % semilogy([-M1:M2],abs(sorted_DNoise3),'bo')
    
   % hold off
    
    
