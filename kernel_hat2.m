

% this function is to compute fourier coefficients of kernel k


function y = kernel_hat2(R,m,k,alpha,X,Y,epsilon)
   
       y = zeros(max(size(X)));

  index  = k^2 - (X+alpha).^2 - (Y*pi/(m*R+epsilon)).^2;

  y(index ~= 0) = (cos(Y(index ~= 0).*pi).*exp(1i*(m*R+epsilon).*sqrt(k^2-(X(index ~= 0)+alpha).^2)) - 1)./...      
                   ((k^2 - (X(index ~= 0)+alpha).^2 - (Y(index ~= 0)*pi/(m*R+epsilon)).^2).*sqrt(2*pi*(2*m*R+2*epsilon)));
  
  y(index == 0) = 1i./(4*abs(Y(index == 0)))*(((m*R+epsilon)/pi)^(1.5));
  
end
