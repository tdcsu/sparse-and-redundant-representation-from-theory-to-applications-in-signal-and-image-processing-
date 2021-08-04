% Figure - 4.2
% =========================================
% This program demonstrates a small part of the proof of the 
% BP performance for the two-ortho case - this code appears in 
% Figure 4.2

n=50; kp=7; kq=9; mu=0.1;
C=[ones(1,kp), ...
    zeros(1,n-kp), ...
    ones(1,kq), ...
    zeros(1,n-kq)];
A=[ones(1,2*n); 
     -ones(1,2*n); 
     eye(n), -ones(n)*mu; 
     -ones(n)*mu, eye(n)]; 
b=[1; -1; zeros(2*n,1)];

x=linprog(-C,A,b,[],[],zeros(2*n,1)); 
plot(x,'*'); 