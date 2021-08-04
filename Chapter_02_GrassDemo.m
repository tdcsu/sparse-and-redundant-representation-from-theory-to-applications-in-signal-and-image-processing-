% Figures 2.3, 2.4, and 2.5
% =========================================
% In this program we create a matrix of size n*m such that its 
% coherence is as small as possible
n=50; m=100; 
D=randn(n,m); 
D=D*diag(1./sqrt(diag(D'*D))); 
[Dout,G,Res]=Chapter_02_DesignFrame(n,m,10000,0.90,0.90,D);

figure(1); clf; 
imagesc([D'*D,ones(m,2),Dout'*Dout; ones(2,2*m+2); 
              ones(m,m+2),abs(Dout'*Dout)]);
colorbar; 
colormap(gray(256))
axis image
axis off;
% print -depsc2 Chapter_02_GrassResults.eps

figure(2); clf
h=semilogx(Res(:,3),'b'); hold on; 
set(h,'LineWidth',2); 
h=semilogx(Res(:,2),'r'); hold on; 
set(h,'LineWidth',2); 
h=semilogx(Res(:,1),'g');
set(h,'LineWidth',2); 
axis([0 10000 0 0.6]); 
legend({'Obtained \mu','mean coherence','Optimal \mu'},1); 
set(gca,'FontSize',14);
% print -depsc2 Chapter_02_GrassConvergence.eps

figure(3); clf; 
gg=D'*D; 
gg=sort((gg(:))); 
gg=gg(1:end-m); 
h=plot(gg); 
set(h,'LineWidth',2); 
hold on; 
gg=Dout'*Dout; 
gg=sort((gg(:))); 
gg=gg(1:end-m); 
h=plot(gg,'r'); 
set(h,'LineWidth',2); 
grid on; 
axis([0 m^2-m -0.6 0.6]);
optmu=sqrt((m-n)/n/(m-1));
h=plot([1,m*(m-1)],[optmu,optmu],'g'); 
set(h,'LineWidth',2); 
legend({'Initial Gram','Final Gram','Optimal \mu'},2); 
h=plot([1,m*(m-1)],[-optmu,-optmu],'g'); 
set(h,'LineWidth',2); 
set(gca,'FontSize',14);
h=plot(gg,'r'); 
set(h,'LineWidth',2); 
% print -depsc2 Chapter_02_GrassResultsVal.eps

