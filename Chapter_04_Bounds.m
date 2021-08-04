% Figure - 4.1
% =========================================
% This program generates the figure that shows the various 
% obtained bounds in the two ortho case

figure(1); clf; 
mu=0.1; 
kp=0.01:0.01:10; 

% The uniqueness reqiurement
kq=1/mu-kp; 
h=plot(kp,kq,'b');
set(h,'LineWidth',2); hold on;

% % The requirement kp*kq<1/mu^2
% kq=1/mu^2./kp
% plot(kp,kq); 
% axis([0 15 0 15]);

% The requirement kp+kq<(sqrt(2)-0.5)/mu
kq=(sqrt(2)-0.5)/mu-kp; 
h=plot(kp,kq,'b-.');
set(h,'LineWidth',2); hold on;

% The exact BP requirement 
kq=(1-mu*kp)./(2*mu^2*kp); 
pos=find(kq<=kp); 
h=plot(kp(pos),kq(pos),'r');
set(h,'LineWidth',2); hold on;

% The OMP requirement
h=plot(0:0.05/mu:0.5/mu,0.5/mu*ones(1,11),'g');
set(h,'LineWidth',2); 
h=plot(0.5/mu*ones(1,11),0:0.05/mu:0.5/mu,'g');
set(h,'LineWidth',2); 

% BP - the second part
h=plot(kq(pos),kp(pos),'r');
set(h,'LineWidth',2); hold on;

% Organizing the graph
grid on;
plot(kp,kp,'k:');
xlabel('k_p');
ylabel('k_q');
%title('The Various Bounds for the Two-Ortho Case');
legend({'Uniqueness bound: 1/\mu','BP bound: 0.91/\mu',...
             'BP exact bound','OMP bound: max(k_p,k_q)=1/2\mu',},1); 
axis image; 
axis([0 10 0 10]); 
% print -depsc2 Chapter_04_BoundsComp2o.eps