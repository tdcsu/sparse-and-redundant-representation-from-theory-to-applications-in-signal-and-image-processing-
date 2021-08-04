% Figure 1.1
% =========================================
% Demonstrating the Lp-Lq relationship and tendency to sparsity

t=[0:1:360]/360*2*pi; 
x=cos(t);
y=sin(t); 

% Case 1
p=2; q=1;

xp=x./(abs(x).^p+abs(y).^p).^(1/p); 
yp=y./(abs(x).^p+abs(y).^p).^(1/p); 
xq=x./(abs(x).^q+abs(y).^q).^(1/q); 
yq=y./(abs(x).^q+abs(y).^q).^(1/q); 

figure(1); clf;
h=plot(xp,yp,'--k'); hold on;
set(h,'LineWidth',2);
h=plot(xq,yq,'k'); 
set(h,'LineWidth',2);
axis image; 
axis([-1.5 1.5 -1.5 1.5]); 
grid on;

factor=sqrt(2);
h=plot(factor*xq,factor*yq,'-.k'); 
set(h,'LineWidth',2);
set(gca,'FontSize',18);
% print -depsc2 Chapter_01_LpLq_demo1.eps

% Case 2
p=1; q=0.5;

xp=x./(abs(x).^p+abs(y).^p).^(1/p); 
yp=y./(abs(x).^p+abs(y).^p).^(1/p); 
xq=x./(abs(x).^q+abs(y).^q).^(1/q); 
yq=y./(abs(x).^q+abs(y).^q).^(1/q); 

figure(1); clf;
h=plot(xp,yp,'--k'); hold on;
set(h,'LineWidth',2);
h=plot(xq,yq,'k'); 
set(h,'LineWidth',2);
axis image; 
axis([-1.5 1.5 -1.5 1.5]); 
grid on;

factor=2;
h=plot(factor*xq,factor*yq,'-.k'); 
set(h,'LineWidth',2);
set(gca,'FontSize',18);
% print -depsc2 Chapter_01_LpLq_demo2.eps


