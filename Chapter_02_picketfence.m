% Figure 2.1
% =========================================
% The PickeFence signal

n=64; 
Pf=zeros(n,1); 
for k=1:sqrt(n):n
    Pf(k)=1; 
end;
figure(1); clf; 
h=plot(1:1:64,Pf,'k'); hold on;
set(h,'LineWidth',2);
h=plot(1:1:64,Pf,'ok'); 
set(h,'MarkerSize',2,'MarkerFaceColor','k');
set(h,'MarkerSize',5);
axis([0 64 -0.1 1.1]);
set(gca,'FontSize',18);
% print -depsc2 Chapter_02_picketfence.eps