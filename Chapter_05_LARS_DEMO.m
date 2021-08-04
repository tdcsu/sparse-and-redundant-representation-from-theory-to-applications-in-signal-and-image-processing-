% Figures - 5.10 and 5.11
% =========================================
% This program presents the LARS solution to the Q problem. 
% We demonstrate several things here:
%       1. Effect of lambda and the optimal value
%       2. Kinds of solutions obtained for various lambda
%       3. The convergence of the IRLS
%       4. The solution path (LARS style)

% Creating the data to test on
load temp.mat % created by running Chapter_15_IRLS

% Alternatively: create the date again differently
% n=100; m=200; S=4; sigma=0.1;
% A=randn(n,m);
% W=sqrt(diag(A'*A));
% for k=1:1:m,
%     A(:,k)=A(:,k)/W(k);
% end;
% x0=zeros(m,1);
% pos=randperm(m);
% x0(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
% b=A*x0+randn(n,1)*sigma;

% The LARS algorithm for varying lambda
Err=zeros(n,1);
Res=zeros(n,1);
[XsolLARS,lambda]=Chapter_05_LARS(A,b,0);
lambda=lambda/2; % missing factor 2 in LASSO formulation
for k=1:1:n
    xout=XsolLARS(k,:)';
    ErrLARS(k)=(xout-x0)'*(xout-x0)/(x0'*x0);
end;

XsolIRLS=zeros(m,n);
ErrIRLS=zeros(100,1);
AA=A'*A;
Ab=A'*b;
for ll=1:1:n
    xout=ones(m,1);
    XX=eye(m);
    for k=1:1:15
        xout=inv(2*lambda(ll)*XX+AA)*Ab;
        XX=diag(1./(abs(xout)+1e-10));
    end;
    Xsol(:,ll)=xout; 
    ErrIRLS(ll)=(xout-x0)'*(xout-x0)/(x0'*x0);
end;

figure(1); clf; 
h=semilogx(lambda,ErrLARS);
set(h,'LineWidth',2); 
hold on;
h=semilogx(lambda,ErrIRLS','r');
set(h,'LineWidth',2); 
h=xlabel('\lambda'); 
set(h,'FontSize',14);
h=ylabel('Normalized L_2-Error');
set(h,'FontSize',14);
set(gca,'FontSize',14);
grid on;
axis([min(lambda) max(lambda) 0 1]); 
legend({'LARS','IRLS'},2);
% print -depsc2 Chapter_05_LARS_lambdaCurve.eps

figure(1); 
axis([0.02    0.2    0.04 0.15]);
% print -depsc2 Chapter_05_LARS_lambdaCurveEnlarged.eps

