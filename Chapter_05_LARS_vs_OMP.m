% Figure - 5.12, 5.13, and 5.14
% =========================================
% This program presents the LARS and OMP approximations 

n=30; m=50; Smax=15; Exper=200; sigma=0.1;

A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 

Er2=zeros(Smax,Exper,3); 
ErS=zeros(Smax,Exper,2); 
Card=zeros(Smax,Exper,2);

for S=1:1:Smax,
    
    h=waitbar(0,'Experimenting ...');
    for experiment=1:1:Exper
        waitbar(experiment/Exper); 
        
        % Generate a test signal of cardinality S
        x0=zeros(m,1);
        pos=randperm(m);
        x0(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
        b=A*x0+randn(n,1)*sigma; 

        % Apply OMP
        xOMP=zeros(m,n); 
        r=b;
        SS=[];
        for k=1:1:n
            Z=abs(A'*r);
            posZ=find(Z==max(Z));
            SS=sort([SS,posZ(1)]);
            r=b-A(:,SS)*pinv(A(:,SS))*b;    
            xOMP(:,k)=zeros(m,1); 
            xOMP(SS,k)=pinv(A(:,SS))*b;
        end;
        choice=find(sqrt(mean((A*xOMP-b*ones(1,n)).^2,1))<=sqrt(2)*sigma,1);
        Sopt=find(xOMP(:,choice))'; 
        ErS(S,experiment,1)=(max(S,length(Sopt))-...
                       length(intersect(Sopt,pos(1:S))))/max(S,length(Sopt)); 
        Er2(S,experiment,1)=mean((xOMP(:,choice)-x0).^2)/(x0'*x0); 
        Card(S,experiment,1)=length(Sopt);
   
        % Apply LARS
        [xLARS,lambda]=Chapter_05_LARS(A,b,0);
        xLARS=xLARS'; 
        choice=find(sqrt(mean((A*xLARS-b*ones(1,n)).^2,1))<=sqrt(2)*sigma,1);
        Sopt=find(xLARS(:,choice))'; 
        ErS(S,experiment,2)=(max(S,length(Sopt))-...
                       length(intersect(Sopt,pos(1:S))))/max(S,length(Sopt)); 
        Er2(S,experiment,2)=mean((xLARS(:,choice)-x0).^2)/(x0'*x0); 
        Card(S,experiment,2)=length(Sopt);
        
        % Projected LARS
        xLARSproj=zeros(m,1); 
        xLARSproj(Sopt)=pinv(A(:,Sopt))*b;
        Er2(S,experiment,3)=mean((xLARSproj-x0).^2)/(x0'*x0); 
    end;
    close(h); 
    
    % displaying the results for the set of tests
    disp([mean(Er2(S,:,1)),mean(Er2(S,:,2)),mean(Er2(S,:,3)),...
             mean(ErS(S,:,1)),mean(ErS(S,:,2))]); 
    
end;

figure(1); clf; 
h=plot(1:1:Smax,mean(Er2(:,:,1),2),'k'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,2),2),'k--'); 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Er2(:,:,3),2),'k:'); 
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Average and Relative L_2-Error'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'OMP','LARS','LARS+Projection'},2);
axis([1 Smax 0 0.02]);
% print -depsc2 Chapter_05_LARSvsOMP_L2.eps

figure(2); clf; 
h=plot(1:1:Smax,mean(ErS(:,:,1),2),'k'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(ErS(:,:,2),2),'k--'); 
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Probability of Error in Support'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'OMP','LARS'},2);
axis([1 15 0 0.6]); 
% print -depsc2 Chapter_05_LARSvsOMP_Supp.eps

figure(3); clf; 
h=plot(1:1:Smax,mean(Card(:,:,1),2),'k'); hold on; 
set(h,'LineWidth',2);
h=plot(1:1:Smax,mean(Card(:,:,2),2),'k--'); 
set(h,'LineWidth',2);
h=plot([0,Smax],[0,Smax],'k:');
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Average Cardinality Estimated'); set(h,'FontSize',14);
set(gca,'FontSize',14);
h=legend({'OMP','LARS','True-cardinality'},2);
axis([1 Smax 1 25]);
% print -depsc2 Chapter_05_LARSvsOMP_Sparse.eps
