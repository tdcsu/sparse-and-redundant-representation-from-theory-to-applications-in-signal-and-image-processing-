% Figure - 9.2
% =========================================
% This program presents the 1D TV analysis operator, and the 
% corresponding Heavyside synthesis one.


% Create the TV matrix
n=100; 
T=zeros(n); 
for k=1:1:99,
    T(k,k)=1; 
    T(k,k+1)=-1; 
end;
T(100,100)=1;

% Compute the equivalent synthesis operator
A=inv(T);

% Show them both
figure(1); clf; 
subplot(1,2,1); 
imagesc(T); axis image; 
colormap(gray(256));
colorbar; 
subplot(1,2,2); 
imagesc(A); axis image; 
colormap(gray(256));
colorbar; 
% print -depsc2 Chapter_09_TV_vs_Heaviside.eps