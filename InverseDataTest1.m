




% Test data

clear all

M1 = 20; M2 = 20; 

k = 1.5;

theta = pi/2;

structure = 1;

q = 0.1;

max = 512;
 
 v = data(max, M1, M2, k, theta, structure, q);
 
for j = 4:-1:1       
    
    NN(j) = 2^(4+j) ;
    
    N = NN(j);    
    
    u = data(N, M1, M2, k, theta, structure, q); 
    
    time(j) = toc;    
            
    relativeError(j) = norm(v-u)/norm(v);

    fprintf('Points: %d \n',N);
    fprintf('Relative error: %f \n',relativeError(j));
    fprintf('Time for numeric computation: %f seconds \n',time(j));
end

save('convergenceTest', 'NN', 'time', 'relativeError')    
loglog(NN,relativeError,'ro');
hold on
loglog(NN,(NN./2).^(-1));

%  figure
%  nearField(1024,60,8*sqrt(2),pi/2,1,0.1);


% nearField2(1024,40,40,8*sqrt(2),pi/2,1,0.1,0);
% 
% figure
% 
% nearField2(1024,40,40,8*sqrt(2),pi/2,1,0.1,40);


