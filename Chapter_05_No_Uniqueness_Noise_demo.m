% Figures - 5.1 and 5.2
% =========================================
% This script demonstrate the lack of uniqueness 
% in solving P_0^{\eps}, first with weak noise, that 
% gives no alternative support

A=[eye(2),[1 -1; 1 1]/sqrt(2)];
x0=[0 0 1 0]';
e=randn(2,1); e=e/norm(e);
b0=A*x0;
eps=0.2; 
b=A*x0+eps*e;

figure(1); clf;
h=plot(b0(1),b0(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 
hold on; 
h=plot(b(1),b(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 

t=0:0.001:2*pi; 
h=plot(b0(1)+eps*cos(t),b0(2)+eps*sin(t),'b');
set(h,'LIneWidth',2); 
axis equal; 
h=plot(b(1)+eps*cos(t),b(2)+eps*sin(t),'r');
set(h,'LIneWidth',2); 

c='rbgc';
C=[255 0 0; 0 0 255; 0 255 0; 0 255 255]/256; 
Outside=[]; 
for k=1:1:4,
    B=A(:,k)*(-10:0.001:10); 
    Flag=sqrt(sum((B-b*ones(1,20001) ).^2,1));
    pos=find(Flag<=eps);
    if ~isempty(pos),
        h=plot(B(1,pos),B(2,pos),'.');
        set(h,'Color',C(k,:),'MarkerSize',15); 
    end;
    pos=find(Flag>eps);
    Outside=[Outside,B(:,pos)];
end;
h=plot(Outside(1,:),Outside(2,:),'.');
set(h,'Color',[0.8 0.8 0.8]); 

axis([-0.2 1.5 -0.2 1.5]); 
legend({'{\bf Ax}_0','{\bf b}',...
            '{v | ||v-{\bf Ax}_0||< \epsilon}',...
            '{v | ||v-{\bf b}||< \epsilon}',...
            '||z \cdot a_3 - b||  < \epsilon','Other sparse solutions'},4); 

h=plot(b0(1),b0(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 
hold on; 
h=plot(b(1),b(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 
% print -depsc2 Chapter_05_NoUniqueness1.eps

% The second part assumes stronger noise, and thus
% alternative supports exist

A=[eye(2),[1 -1; 1 1]/sqrt(2)];
x0=[0 0 1 0]';
e=randn(2,1); e=e/norm(e);
b0=A*x0;
eps=0.6; 
b=A*x0+eps*e;

figure(1); clf;
h=plot(b0(1),b0(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 
hold on; 
h=plot(b(1),b(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 

t=0:0.001:2*pi; 
h=plot(b0(1)+eps*cos(t),b0(2)+eps*sin(t),'b');
set(h,'LIneWidth',2); 
axis equal; 
h=plot(b(1)+eps*cos(t),b(2)+eps*sin(t),'r');
set(h,'LIneWidth',2); 

c='rbgc';
C=[255 0 0; 0 0 255; 0 255 0; 0 255 255]/256; 
Outside=[]; 
for k=1:1:4,
    B=A(:,k)*(-10:0.001:10); 
    Flag=sqrt(sum((B-b*ones(1,20001) ).^2,1));
    pos=find(Flag<=eps);
    if ~isempty(pos),
        h=plot(B(1,pos),B(2,pos),'.');
        set(h,'Color',C(k,:),'MarkerSize',15); 
    end;
    pos=find(Flag>eps);
    Outside=[Outside,B(:,pos)];
end;
h=plot(Outside(1,:),Outside(2,:),'.');
set(h,'Color',[0.8 0.8 0.8]); 

axis([-0.2 1.5 -0.2 1.5]); 
legend({'{\bf Ax}_0','{\bf b}',...
            '{v | ||v-{\bf Ax}_0||< \epsilon}',...
            '{v | ||v-{\bf b}||< \epsilon}',...
            '||z \cdot a_1 - b||  < \epsilon',...
            '||z \cdot a_2 - b||  < \epsilon',...
            '||z \cdot a_3 - b||  < \epsilon',...
            '||z \cdot a_4 - b||  < \epsilon',...
            'Other sparse solutions'},4); 
        
h=plot(b0(1),b0(2),'o'); 
set(h,'MarkerSize',10,'MarkerFaceColor','b'); 
hold on; 
h=plot(b(1),b(2),'sr'); 
set(h,'MarkerSize',10,'MarkerFaceColor','r'); 
% print -depsc2 Chapter_05_NoUniqueness2.eps

