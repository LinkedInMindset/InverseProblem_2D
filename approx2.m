



function coeff = approx2(epsilon)

% cut-off function is approximated by polynomial of degree 6 in two period
% (-2*rho,2*rho) and (2*rho,6*rho)
%rho = 1;




A = [-(epsilon/2)^7 +(epsilon/2)^6 -(epsilon/2)^5 (epsilon/2)^4 -(epsilon/2)^3 (epsilon/2)^2 -(epsilon/2) 1; 
   7*(epsilon/2)^6 -6*(epsilon/2)^5 5*(epsilon/2)^4 -4*(epsilon/2)^3 3*(epsilon/2)^2 -2*(epsilon/2) 1 0;-42*(epsilon/2)^5 30*(epsilon/2)^4 -20*(epsilon/2)^3 12*(epsilon/2)^2 -6*(epsilon/2) 2 0 0;
   210*(epsilon/2)^4 -120*(epsilon/2)^3 60*(epsilon/2)^2 -24*(epsilon/2) 6 0 0 0;(epsilon/2)^7 (epsilon/2)^6 (epsilon/2)^5 (epsilon/2)^4 (epsilon/2)^3 (epsilon/2)^2 (epsilon/2) 1;
  7*(epsilon/2)^6 6*(epsilon/2)^5 5*(epsilon/2)^4 4*(epsilon/2)^3 3*(epsilon/2)^2 2*(epsilon/2) 1 0;
   42*(epsilon/2)^5 30*(epsilon/2)^4 20*(epsilon/2)^3 12*(epsilon/2)^2 6*(epsilon/2) 2 0 0; 210*(epsilon/2)^4 120*(epsilon/2)^3 60*(epsilon/2)^2 24*(epsilon/2) 6 0 0 0];

b = [1; 0; 0; 0; 0; 0; 0; 0];

coeff = A\b;














%  m = 4; rho = 1;
% 
%  x=(0:1/1000:m*rho + epsilon);
%  
%  z=zeros(1,length(x));
% 
%  z(x<=m*rho) = 1;
%  
%  z(x>=m*rho) = coeff(1)*(x(x>=m*rho)-m*rho - epsilon/2).^7 +coeff(2)*(x(x>=m*rho)-m*rho - epsilon/2).^6 + coeff(3)*(x(x>=m*rho)-m*rho - epsilon/2).^5 + coeff(4)*(x(x>=m*rho)-m*rho - epsilon/2).^4 ...
%         +coeff(5)*(x(x>=m*rho)-m*rho - epsilon/2).^3 + coeff(6)*(x(x>=m*rho)-m*rho - epsilon/2).^2 + coeff(7)*(x(x>=m*rho)-m*rho - epsilon/2) + coeff(8);
%    
%   plot(x,z)


%  x=[-2*rho:1/100000:6*rho];
%  
%   z=zeros(1,length(x));
%  
%  
%   z(x<=0) = coeff(1)*(-x(x<=0)-3*rho/2).^8 +coeff(2)*(-x(x<=0)-3*rho/2).^7 +coeff(3)*(-x(x<=0)-3*rho/2).^6 + coeff(4)*(-x(x<=0)-3*rho/2).^5 ... 
%       +coeff(5)*(-x(x<=0)-3*rho/2).^4 + coeff(6)*(-x(x<=0)-3*rho/2).^3 + coeff(7)*(-x(x<=0)-3*rho/2).^2 + coeff(8)*(-x(x<=0)-3*rho/2) + 1;
%   
%   z(x>=0 & x<=2*rho) = coeff(1)*(x(x>=0 & x<=2*rho)-3*rho/2).^8 +coeff(2)*(x(x>=0 & x<=2*rho)-3*rho/2).^7 + coeff(3)*(x(x>=0 & x<=2*rho)-3*rho/2).^6 + coeff(4)*(x(x>=0 & x<=2*rho)-3*rho/2).^5 ...
%        + coeff(5)*(x(x>=0 & x<=2*rho)-3*rho/2).^4 + coeff(6)*(x(x>=0 & x<=2*rho)-3*rho/2).^3 + coeff(7)*(x(x>=0 & x<=2*rho)-3*rho/2).^2 + coeff(8)*(x(x>=0 & x<=2*rho)-3*rho/2) + 1;
%   
%   z(x>2*rho & x<=4*rho) = coeff(1)*(-x(x<0)-3*rho/2).^8 + coeff(2)*(-x(x<0)-3*rho/2).^7 + coeff(3)*(-x(x<0)-3*rho/2).^6 + coeff(4)*(-x(x<0)-3*rho/2).^5 ...
%       + coeff(5)*(-x(x<0)-3*rho/2).^4 + coeff(6)*(-x(x<0)-3*rho/2).^3 + coeff(7)*(-x(x<0)-3*rho/2).^2 + coeff(8)*(-x(x<0)-3*rho/2) + 1;
%   
%   z(x>=4*rho) = coeff(1)*(x(x>=0 & x<=2*rho)-3*rho/2).^8 +coeff(2)*(x(x>=0 & x<=2*rho)-3*rho/2).^7 + coeff(3)*(x(x>=0 & x<=2*rho)-3*rho/2).^6 + coeff(4)*(x(x>=0 & x<=2*rho)-3*rho/2).^5 ...
%       +coeff(5)*(x(x>=0 & x<=2*rho)-3*rho/2).^4 + coeff(6)*(x(x>=0 & x<=2*rho)-3*rho/2).^3 + coeff(7)*(x(x>=0 & x<=2*rho)-3*rho/2).^2 + coeff(8)*(x(x>=0 & x<=2*rho)-3*rho/2) + 1;
%   
%   plot(x,z)