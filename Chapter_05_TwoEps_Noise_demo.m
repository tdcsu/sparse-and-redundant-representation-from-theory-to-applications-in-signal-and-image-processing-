% Figure - 5.3
% =========================================
% This script demonstrate the distance of 2\eps between Ax1 and Ax2. 

eps = 0.2; 
b = [1; 2]; 
e = randn(2,1); 
e = e/norm(e);
Ax1 = b+eps*e; 
e = randn(2,1); 
e = e/norm(e);
Ax2 = b+eps*e; 

figure(1); clf;
h=plot(b(1),b(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 
hold on; 
h=plot(Ax1(1),Ax1(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 
h=plot(Ax2(1),Ax2(2),'sg'); 
set(h,'MarkerSize',10,'MarkerFaceColor','g'); 

t=0:0.001:2*pi; 
h=plot(b(1)+eps*cos(t),b(2)+eps*sin(t),'b');
set(h,'LIneWidth',2); 

axis equal; 
axis([0.7 1.3 1.7 2.3]); 
axis off; 

disp('Mark the point where the radius should touch the circle');
[xx,yy]=ginput(1);
h=quiver(b(1),b(2),xx-b(1),yy-b(2),1);
set(h,'LineWidth',2,'MaxHeadSize',10,'Color',[0 0 0]);
h=plot(b(1),b(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 

h=quiver(Ax1(1),Ax1(2),Ax2(1)-Ax1(1),Ax2(2)-Ax1(2),0.95);
set(h,'LineWidth',2,'MaxHeadSize',10,'Color',[0 0 0]);
h=quiver(Ax2(1),Ax2(2),Ax1(1)-Ax2(1),Ax1(2)-Ax2(2),0.95);
set(h,'LineWidth',2,'MaxHeadSize',10,'Color',[0 0 0]);
h=plot(Ax1(1),Ax1(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 
h=plot(Ax2(1),Ax2(2),'sg'); 
set(h,'MarkerSize',10,'MarkerFaceColor','g'); 


disp('Position the text b in the center');
h=gtext('b'); set(h,'FontSize',14); 
disp('Position the text Ax_1 near the red box');
h=gtext('Ax_1'); set(h,'FontSize',14); 
disp('Position the text Ax_2 near the green box');
h=gtext('Ax_2'); set(h,'FontSize',14); 
disp('Position the text \epsilon on the radius arrow');
h=gtext('\epsilon'); set(h,'FontSize',14); 
disp('Position the text ||Ax_1 - Ax_2||_2 near the two-sided-arrow');
h=gtext('||Ax_1 - Ax_2||_2'); set(h,'FontSize',14); 

% print -depsc2 Chapter_05_NoisyDistance.eps