
function  [picture,q11_] = inverseSolver2TikhonovWN_FullSpace(N, M1, M2, k, m , R, R_meas, theta, accuracy, q11, q22, epsilon, structure, lastEigenvalues, noiselevel, gamma)


      h1 = 2*pi/N; 
      h2 = (2*m*R + 2*epsilon)/N;   
   alpha = k*cos(theta);    
    
   % Compute the coefficients for implementing VFFT and IVFFT
   B1=ones(N,1)*exp(-2*pi*1i*(1-N/2)*(1:N)/(N));
   B2=ones(N,1)*exp(-2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));
   C1=ones(N,1)*exp(2*pi*1i*(1-N/2)*(1:N)/(N));
   C2=ones(N,1)*exp(2*pi*1i*(((1-N/2)*(1:N)+(N)^2/4-1))/(N));    

   
  [XN,YN] = meshgrid(-N/2+1:N/2); % X-left-right, Y-down-up  
  
  % Compute the contrast q numerically
  q11_ = zeros(N); q22_ = zeros(N);
  
  if strcmp(structure,'sin')==1      
       index_sin = (YN*h2 < -cos(XN*h1*2)/4 +0.25) & (YN*h2 > -0.5);
       q11_(index_sin) = q11; 
       q22_(index_sin) = q22;
       %q22_(abs(YN*h2) < -cos(XN*h1)/2 +1) = q22;  
       
  elseif strcmp(structure,'sin_smooth')==1      
       cutoff = zeros(N);
       index_cut = (abs(YN*h2) < 0.5);
       cutoff(index_cut) = exp(1 - 1./(1 - ((YN(index_cut)*h2).^2)/0.5^2));
       index_sin = (YN*h2 < -cos(XN*h1*2)/4+0.25) & (YN*h2 > -0.5);
       q11_(index_sin) = q11*cutoff(index_sin); 
       q22_(index_sin) = q22*cutoff(index_sin);    
       
  elseif strcmp(structure,'circle')==1
       r = 0.5;
       index_cir = ((YN*h2).^2 + (XN*h1).^2 <= r^2);
       q11_(index_cir) = q11*exp(1 - 1./(1 - ((XN(index_cir)*h1).^2+(YN(index_cir)*h2).^2)/r^2));   
       q22_ = 0.5*q11_;
       
  elseif strcmp(structure,'circle2')==1
       r1 = 0.5; x1=-pi/2;y1=0;x2=pi/2;y2=0;  r2 = 0.5;
       index_cir1 = ((XN*h1-x1).^2+(YN*h2-y1).^2 <= r1^2);
       index_cir2 = ((XN*h1-x2).^2+(YN*h2-y2).^2 <= r2^2);
       q11_(index_cir1) = q11*exp(1 - 1./(1 - ((XN(index_cir1)*h1-x1).^2+(YN(index_cir1)*h2-y1).^2)/r1^2));  
       q11_(index_cir2) = q11*exp(1 - 1./(1 - ((XN(index_cir2)*h1-x2).^2+(YN(index_cir2)*h2-y2).^2)/r2^2));   

       q22_ = 0.5*q11_;
       
  elseif strcmp(structure,'circle3')==1
       r1 = 0.5; x1=-pi-1+2*pi/3;y1=0;
       r2 = 0.5; x2=-pi-1+4*pi/3;y2=0;  
       r3 = 0.5; x3=-pi-1+6*pi/3;y3=0;
       index_cir1 = ((XN*h1-x1).^2+(YN*h2-y1).^2 <= r1^2);
       index_cir2 = ((XN*h1-x2).^2+(YN*h2-y2).^2 <= r2^2);
       index_cir3 = ((XN*h1-x3).^2+(YN*h2-y3).^2 <= r3^2);

       q11_(index_cir1) = q11*exp(1 - 1./(1 - ((XN(index_cir1)*h1-x1).^2+(YN(index_cir1)*h2-y1).^2)/r1^2));  
       q11_(index_cir2) = q11*exp(1 - 1./(1 - ((XN(index_cir2)*h1-x2).^2+(YN(index_cir2)*h2-y2).^2)/r2^2));   
       q11_(index_cir3) = q11*exp(1 - 1./(1 - ((XN(index_cir3)*h1-x3).^2+(YN(index_cir3)*h2-y3).^2)/r3^2));   

       q22_ = 0.5*q11_;     
       
  elseif strcmp(structure,'cube_ball')==1        
        r1 = 0.1; r = 0.4;
       
       ball = (XN*h1).^2 + (YN*h2-0.4).^2 <= r^2;
       q11_(ball) = q11*exp(1 - 1./(1 - ((XN(ball)*h1).^2+(YN(ball)*h2-0.4).^2)/r^2));
       f1 = zeros(N); 
       f1(abs(YN*h2)<r1) = exp(1-r1^2./(r1^2 - (YN(abs(YN*h2)<r1)*h2).^2));
       q11_(abs(YN*h2)<r1) = q11*f1(abs(YN*h2)<r1);
       q22_ = 0.5*q11_;   
       
  elseif strcmp(structure,'ellip')==1        
       q11_((YN*h2).^2/(0.75^2) + (XN*h1).^2/(2.5^2) <= 1) = q11;   
       q22_ = 0.5*q11_; 
       
  elseif strcmp(structure,'halfEllip')==1        
       q11_(((YN*h2).^2/(0.75^2) + (XN*h1+3).^2/(6^2) <= 1)& XN*h1 >=-3) = q11;
       q22_ = 0.5*q11_;
       
  elseif strcmp(structure,'cube_smooth')==1        
        r1 = 2.5; r2 = 0.35;
       f1 = zeros(N); f2 = zeros(N);
       f1(abs(XN*h1)<r1) = exp(1-r1^2./(r1^2 - (XN(abs(XN*h1)<r1)*h1).^2));
       f2(abs(YN*h2)<r2) = exp(1-r2^2./(r2^2 - (YN(abs(YN*h2)<r2)*h2).^2));
       q11_ = q11*f1.*f2;
       q22_ = 0.5*q11_;
       
  elseif strcmp(structure,'cube_smooth2')==1        
        r1 = 1.8; r2 = 0.25;r3 = 0.25;
       f1 = zeros(N); f2 = zeros(N); f3 = zeros(N);
       f1(abs(XN*h1)<r1) = exp(1-r1^2./(r1^2 - (XN(abs(XN*h1)<r1)*h1).^2));
       f2(abs(YN*h2-0.1)<r2) = exp(1-r2^2./(r2^2 - (YN(abs(YN*h2-0.1)<r2)*h2-0.1).^2));
       f3(abs(YN*h2+0.25)<r3) = exp(1-r3^2./(r3^2 - (YN(abs(YN*h2+0.25)<r3)*h2+0.25).^2));
       
       q11_ = q11*(f1.*f2 + f3)/2;
       
       q22_ = 0.5*q11_;

  elseif strcmp(structure,'cube1')==1        
       q11_(((XN*h1>=-pi & XN*h1<=-pi/2) | (XN*h1<=pi & XN*h1>=pi/2)) & abs(YN*h2)<=1.5) = 0.5;
       q11_((XN*h1>-pi/2 & XN*h1<pi/2) & abs(YN*h2)<=0.5) = 1-3*1i;   
       q22_ = q11_;
       
  elseif strcmp(structure,'cube2')==1 

       q11_(((XN*h1<-3*pi/4)|(XN*h1>3*pi/4))  & (YN*h2 <=1 & YN*h2 >=-0.5)) = q11;
       q11_((XN*h1>=-pi/2 & XN*h1<=pi/2) & (YN*h2 <=0 & YN*h2>=-0.5)) = q11;        
       q11_((XN*h1>=pi/2 & XN*h1<=3*pi/4) & (YN*h2 <=4/pi*XN*h1-2 & YN*h2 >=-0.5)) = q11; 
       q11_((XN*h1>=-3*pi/4 & XN*h1<=-pi/2) & (YN*h2 <=-4/pi*XN*h1-2 & YN*h2>=-0.5)) = q11; 
       
       q22_ = q11_;
  elseif strcmp(structure,'cube3')==1        
       q11_(((XN*h1>=-pi & XN*h1<=-pi/2) | (XN*h1<=pi & XN*h1>=pi/2)) & ((YN*h2<=1.5) & (YN*h2>=-0.5))) = 1;
       q11_((XN*h1>-pi/2 & XN*h1<pi/2) & ((YN*h2<=0.5) & (YN*h2>=-1.5))) = 1;   
       q22_ = q11_;
  elseif strcmp(structure,'cube4')==1        
       q11_(abs(XN*h1)<=2.5 & YN*h2 <=0.5 & YN*h2 >=-0.5) = 0.5;       
       q22_ = q11_;
  elseif strcmp(structure,'cube5')==1        
       q11_(((XN*h1>=-pi & XN*h1<=-pi/2) | (XN*h1<=pi & XN*h1>=pi/2)) & ((YN*h2<=1) & (YN*h2>=-0.5))) = 0.5;
       q11_((XN*h1>-pi/2 & XN*h1<pi/2) & (YN*h2<=0) & (YN*h2>=-0.5)) = 0.5;   
       q22_ = q11_;
  elseif strcmp(structure,'cube6')==1        
       q11_(((abs(XN*h1+pi/2)<=0.5) | (abs(XN*h1-pi/2)<=0.5)) & ((YN*h2<=0.5) & (YN*h2>=0))) = 0.4;
       q11_( (YN*h2<=0) & (YN*h2>=-0.5)) = 0.25;   
       q11_ = smooth2a(q11_,5,5);
       q22_ = 0.5*q11_;     
  elseif  strcmp(structure,'cube7')==1      
       a1=-pi-1+2*pi/3;
       a2=-pi-1+4*pi/3;  
       a3=-pi-1+6*pi/3;
       q11_(((abs(XN*h1-a1)<=0.5) | (abs(XN*h1-a2)<=0.5) | (abs(XN*h1-a3)<=0.5)) & ((YN*h2<=0.5) & (YN*h2>=-0.5))) = 0.5;
       q11_( (YN*h2<=0) & (YN*h2>=-0.5)) = 0.5;   
       q22_ = 0.5*q11_; 
       
   elseif  strcmp(structure,'cube8')==1      
       a1=-pi-1+2*pi/3;
       a2=-pi-1+4*pi/3;  
       a3=-pi-1+6*pi/3;
       r1 = 0.5; r2 = 0.5;
       q11_(((abs(XN*h1-a1)<=r1) | (abs(XN*h1-a2)<=r1) | (abs(XN*h1-a3)<=r1)) & ((YN*h2<=r2) & (YN*h2>=-r2))) = 0.5;
       q22_ = 0.5*q11_;    
  end    
  
  
%  Compute WN from direct solver

%  noise = unifrnd(-1,1,[2*(M1+M2+1),2*(M1+M2+1)]) + 1i*unifrnd(-1,1,[2*(M1+M2+1),2*(M1+M2+1)]);
  
    noise = rand(2*(M1+M2+1)) + 1i*rand(2*(M1+M2+1));


%   noise = rand(2*(M1+M2+1));  


 if noiselevel == 0

     [WN,w_starPlus,w_starMinus] = data_CollocationScheme_2D(N, M1, M2, k, m , R, R_meas, theta, accuracy, q11_, q22_, epsilon, B1, B2, C1, C2);  
   
   ReWN = (WN + conj(permute(WN,[1 2 4 3])))/2;
   ImWN = (WN - conj(permute(WN,[1 2 4 3])))/(2*1i);    
   
      
  % chuyen doi dang cua mang ReWN de tim tri rieng va vector rieng 
   
  ImWN1_reshape = reshape(ImWN(:,:,1,:),[(M1+M2+1) 2*(M1+M2+1)]);
  ImWN2_reshape = reshape(ImWN(:,:,2,:),[(M1+M2+1) 2*(M1+M2+1)]);
  
   ImWN_reshape = cat(1,ImWN1_reshape,ImWN2_reshape);  
   
     
  ReWN1_reshape = reshape(ReWN(:,:,1,:),[(M1+M2+1) 2*(M1+M2+1)]);
  ReWN2_reshape = reshape(ReWN(:,:,2,:),[(M1+M2+1) 2*(M1+M2+1)]);

  
  ReWN_reshape = cat(1,ReWN1_reshape,ReWN2_reshape);

%   ReWN_reshape = reshape(ReWN,[2*(M1+M2+1) 2*(M1+M2+1)]);
%   ImWN_reshape = reshape(ImWN,[2*(M1+M2+1) 2*(M1+M2+1)]);
  
  [V1,D1] = eig(ReWN_reshape);  
 
   
  % gia tri tuyet doi cua ma tran duoc dinh nghia thong qua tri rieng va vector rieng 
    WN_reshape = ReWN_reshape + 1i*ImWN_reshape;
    absReWN_reshape = V1*abs(D1)*inv(V1);
  
    
  % tinh WNshape duoi dang mang 6 chieu
  % WNshape = reshape(absReWN_reshape,[M1+M2+1 M1+M2+1 M1+M2+1 M1+M2+1 4 4]) + ImWN;      
  
  % WNshape = absReWN_reshape_reshape + ImWN;  
  WNshape = absReWN_reshape - ImWN_reshape;         
  [V2,D2] = eig(WNshape);  
  
  squareRoot_WNshape = V2*sqrt(real(D2))*inv(V2);
  
   noise = rand(2*(M1+M2+1));  
   
   squareRoot_WNshapeNoise = squareRoot_WNshape + noiselevel*noise*norm(squareRoot_WNshape)/norm(noise) ;
   
   
      if strcmp(structure,'circle2')==1    
       save('circle2.mat');
       elseif strcmp(structure,'sin')==1 
       save('sin.mat');
       elseif strcmp(structure,'circle3')==1 
       save('circle3.mat');
       elseif strcmp(structure,'cube8')==1 
       save('cube8.mat');
       elseif strcmp(structure,'cube_ball')==1 
       save('cube_ball.mat');
       end
  
 else
       if strcmp(structure,'circle2')==1    
       load('circle2.mat','WN_reshape','squareRoot_WNshape','w_starPlus','w_starMinus');
       elseif strcmp(structure,'circle3')==1 
       load('circle3.mat','WN_reshape','squareRoot_WNshape','w_starPlus','w_starMinus');
       elseif strcmp(structure,'sin')==1 
       load('sin.mat','WN_reshape','squareRoot_WNshape','w_starPlus','w_starMinus');
       elseif strcmp(structure,'cube8')==1 
       load('cube8.mat','WN_reshape','squareRoot_WNshape','w_starPlus','w_starMinus');
       elseif strcmp(structure,'cube_ball')==1 
       load('cube_ball.mat','WN_reshape','squareRoot_WNshape','w_starPlus','w_starMinus');
       end
       
    % DAY LA TEST VOI NOISE THEM VAO WN DE TANG TINH THUC TE CUA NOISE-TEST
   % CAC KY HIEU VAN DE NGUYEN NHU LUC THEM NOISE VAO squareroot_WNsharpe NHU O PHIA TREN     
   squareRoot_WNshapeNoise = WN_reshape + noiselevel*noise*norm(WN_reshape)/norm(noise) ;    
 end
   
   
   
   
  %[VNoise3 DNoise3] = eig(squareRoot_WNshapeNoise*squareRoot_WNshapeNoise');

  [UNoise3,DNoise3,V] = svd(squareRoot_WNshapeNoise);
  
  UNoise3 = reshape(UNoise3,[M1+M2+1 1 2 2*(M1+M2+1)]);
  DNoise3 = sum(DNoise3);
  
    
  index = YN(YN*h2<=R + 1 & YN*h2>=-R-1) + N/2 + N*(XN(YN*h2<=R + 1 & YN*h2>=-R-1) + N/2 - 1);  
  %indice of tested domain
  
  size(index);
      
  %index1 = find(k^2 > (alpha + (-M1:M2)).^2);      
  %index1 - M1 - 1  % indice of propagating modes  
  
  
  [sorted_DNoise3,index3] = sort(abs(DNoise3),'descend'); % sort the eigenvalues in descending order  
  %V2 = V2(:,index3); % eigenvectors in corresponding order  
  
  picture = ones(N);  
  
          
  alphaj = alpha + (-M1:M2);
  
    betaj = sqrt(k^2-alphaj.^2); 
    
    
   find(k^2 > alphaj.^2)  % indice of propagating modes  
   
   index2 = (1 : 1 : 2*(M1+M2+1) - lastEigenvalues);
         
%    if noiselevel~= 0
%        
%       gamma = 0.2;
%    else
%       gamma = 0;
%    end
   
     for t = 1:length(index)   
         
          rj_plus = 1i./(4*pi.*betaj).*exp(-1i*(alphaj.*XN(index(t))*h1 + betaj.*(YN(index(t))*h2-R_meas)));        
          rj_minus = 1i./(4*pi.*betaj).*exp(-1i*(alphaj.*XN(index(t))*h1 - betaj.*(YN(index(t))*h2+R_meas))); 
          
                          
          Wr_1 = 4*pi*w_starPlus.*(rj_plus + rj_minus);
          Wr_2 = 4*pi*w_starMinus.*(rj_plus - rj_minus);
                                                              
          
                product1 = conj(UNoise3(:,:,1,:)).*repmat(Wr_1.',[1 1 1 2*(M1+M2+1)]);
                product2 = conj(UNoise3(:,:,2,:)).*repmat(Wr_2.',[1 1 1 2*(M1+M2+1)]);
                          
%                 product1 = conj(UNoise3(:,:,1,:)).*repmat(rj_plus.',[1 1 1 2*(M1+M2+1)]);
%                 product2 = conj(UNoise3(:,:,2,:)).*repmat(rj_minus.',[1 1 1 2*(M1+M2+1)]);   
          
            innerProduct = abs(sum(product1 + product2)).^2;          %((sum(DNoise3(:,index2)) + gamma).^2./sum(DNoise3(:,index2)))
               numerator = cat(1,innerProduct(:)).';
               numerator = numerator(index3);
                
            %   gamma = fzero(@(x) sum(numerator./((DNoise3.^2 + x).^2./(x^2 - (noiselevel)^2*DNoise3.^2))), 0);              

               picture(index(t)) = sum(numerator(:,index2)./((sorted_DNoise3(:,index2).^2 + gamma).^2./sorted_DNoise3(:,index2).^2));
          
 %   picture(index(t)) = sum(numerator(:,index2)./sorted_DNoise3(:,index2).^2);          
                    
     end    
                
    
     