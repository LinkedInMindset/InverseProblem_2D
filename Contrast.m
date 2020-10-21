

function q_N = Contrast(N,X,Y,q,R,epsilon)

%R=1; N=256; epsilon =0.4; q=0.5;

%[X Y] =meshgrid(-N/2+1:N/2);


  h1 = 2*pi/N; h2 = (4*R+2*epsilon)/N;

  q_N = zeros(N);  
 
 
  q_N(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<=0 & Y*h2>=-1/2) = q;   
 
  
  q_N(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<=-1/2 & Y*h2>=-1) = 2*q*Y(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<-1/2 & Y*h2>-1)*h2 + 2*q   ;
  
  
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=0 & Y*h2>=-1/2) = 4*q*X(X*h1>=0 & X*h1<=pi/4 & Y*h2<=0 & Y*h2>=-1/2)*h1/pi  ; 
  
  
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=2*X*h1/pi-1 & Y*h2>=-1) = 2*q*Y(X*h1>=0 & X*h1<=pi/4 & Y*h2<=2*X*h1/pi-1 & Y*h2>=-1)*h2 + 2*q   ;
  
  
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=-1/2 & Y*h2>=2*X*h1/pi-1) = 4*q*X(X*h1>=0 & X*h1<=pi/4 & Y*h2<=-1/2 & Y*h2>=2*X*h1/pi-1)*h1/pi   ;
 
  
  
  
  q_N(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=0 & Y*h2>=-1/2) = q;   
 
  
  q_N(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=-1/2 & Y*h2>=-1) = 2*q*Y(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=-1/2 & Y*h2>=-1)*h2 + 2*q   ;
  
  
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=0 & Y*h2>=-1/2) = 4*q*(pi-X(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=0 & Y*h2>=-1/2)*h1)/pi  ; 
  
   
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=2*(pi-X*h1)/pi-1 & Y*h2>=-1) = 2*q*Y(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=2*(pi-X*h1)/pi-1 & Y*h2>=-1)*h2 + 2*q   ;
   
   
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=-1/2 & Y*h2>=2*(pi-X*h1)/pi-1) = 4*q*(pi-X(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=-1/2 & Y*h2>=2*(pi-X*h1)/pi-1)*h1)/pi   ;
  
   
  
  
  
  q_N(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=1/2 & Y*h2>=0) = q;   
 
  
  q_N(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=1 & Y*h2>=1/2) = 2*q*(R-Y(X*h1>=pi/2 & X*h1<=3*pi/4 & Y*h2<=1 & Y*h2>=1/2)*h2)   ;
  
  
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=1/2 & Y*h2>=0) = 4*q*(pi-X(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=1/2 & Y*h2>=0)*h1)/pi  ; 
   
    
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=2*X*h1/pi-1 & Y*h2>=R/2) = 4*q*(pi-X(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=2*X*h1/pi-1 & Y*h2>=R/2)*h1)/pi  ;
   
    
  q_N(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=1 & Y*h2>=2*X*h1/pi-1) =  2*q*(R-Y(X*h1>=3*pi/4 & X*h1<=pi & Y*h2<=1 & Y*h2>=2*X*h1/pi-1)*h2)   ;
   
  
  
   
  q_N(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<=1/2 & Y*h2>=0) = q;   
 
  
  q_N(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<=1 & Y*h2>=1/2) = 2*q*(R-Y(X*h1>=pi/4 & X*h1<=pi/2 & Y*h2<=1 & Y*h2>=1/2)*h2)  ;
   
  
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=R/2 & Y*h2>=0) = 4*q*X(X*h1>=0 & X*h1<=pi/4 & Y*h2<=R/2 & Y*h2>=0)*h1/pi ; 
   
  
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=2*(pi-X*h1)/pi-1 & Y*h2>=1/2) = 4*q*X(X*h1>=0 & X*h1<=pi/4 & Y*h2<=2*(pi-X*h1)/pi-1 & Y*h2>=1/2)*h1/pi ;  
   
   
  q_N(X*h1>=0 & X*h1<=pi/4 & Y*h2<=1 & Y*h2>=2*(pi-X*h1)/pi-1) = 2*q*(R-Y(X*h1>=0 & X*h1<=pi/4 & Y*h2<=1 & Y*h2>=2*(pi-X*h1)/pi-1)*h2)   ;
  
  
   
  q_N(X*h1>-pi & X*h1<=0 & Y*h2<=1 & Y*h2>=-1) = q_N(X*h1>0 & X*h1<=pi & Y*h2<=1 & Y*h2>=-1) ;
    
   
 
%      mesh((-N/2+1:N/2)*h1,Y(abs(Y(:,1)*h2)<=3)*h2,real(q_N(Y(abs(Y(:,1)*h2)<=3)+N/2,1:N))); % plot the contrast
%      view(2); axis tight
%      axis square 
%      colorbar
     
     
     
     