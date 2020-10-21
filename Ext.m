
% N < M


function y = Ext(x,N,M)

y = zeros(M);

[X Y] = meshgrid(-N/2+1:N/2);

y(Y + M/2 + M*(X-1 + M/2)) = x;