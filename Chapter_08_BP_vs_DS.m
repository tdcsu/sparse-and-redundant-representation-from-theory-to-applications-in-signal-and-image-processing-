% Figures - 8.1 and 8.2
% =========================================
% This program compares the Dantzig Selector to the Basis Pursuit

function []=Chapter_08_BP_vs_DS()

% This program presents the LARS and OMP approximations 

n=50; m=80; S=10; Exper=200; sigma=0.05;
A=dctmtx(m);
pos=randperm(m);
A=A(sort(pos(1:n)),:);
% A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 
AA=A'*A; 
    
Best_BP=0;
Best_DS=0; 

Er2_LARS=zeros(n,Exper);
Er2_LARS_proj=zeros(n,Exper);
ErA2_LARS=zeros(n,Exper);
ErA2_LARS_proj=zeros(n,Exper);
Card_LARS=zeros(n,Exper);

Er2_DSc=cell(n,Exper);
Er2_DS_projc=cell(n,Exper);
ErA2_DSc=cell(n,Exper);
ErA2_DS_projc=cell(n,Exper);

Er2_DS=zeros(n,Exper);
Er2_DS_proj=zeros(n,Exper);
ErA2_DS=zeros(n,Exper);
ErA2_DS_proj=zeros(n,Exper);
Card_DS=zeros(n,Exper);

Er2_ORACLE=zeros(Exper,1);
ErA2_ORACLE=zeros(Exper,1);

for experiment=1:1:Exper
    disp(['XXXXXXXXXXX Experiment number ====>',num2str(experiment)]);
    
    % Generate a test signal of cardinality S
    x0=zeros(m,1);
    pos=randperm(m);
    Sin=pos(1:S);
    x0(pos(1:S))=sign(randn(S,1)).*randn(S,1);
    b=A*x0+randn(n,1)*sigma;
    
    % Apply the oracle
    xOR=zeros(m,1);
    xOR(Sin)=pinv(A(:,Sin))*b;
    Er2_ORACLE(experiment)=mean((xOR-x0).^2)/mean(x0.^2);
    ErA2_ORACLE(experiment)=mean((A*xOR-A*x0).^2)/mean((A*x0).^2);    
    
    % Apply BP (by LARS)
    [xLARS,lambdaLARS]=Chapter_05_LARS(A,b,0); 
    xLARS=xLARS';
    for k=1:1:n
        Sout=find(abs(xLARS(:,k)));
        Er2_LARS(k,experiment)=mean((xLARS(:,k)-x0).^2)/mean(x0.^2);
        ErA2_LARS(k,experiment)=mean((A*xLARS(:,k)-A*x0).^2)/mean((A*x0).^2);
        xLARS_proj=zeros(m,1); 
        xLARS_proj(Sout)=pinv(A(:,Sout))*b;
        Er2_LARS_proj(k,experiment)=mean((xLARS_proj-x0).^2)/mean(x0.^2); 
        ErA2_LARS_proj(k,experiment)=mean((A*xLARS_proj-A*x0).^2)/mean((A*x0).^2); 
        Card_LARS(k,experiment)=length(Sout);
    end;
    
    % Apply DS (using Matlab's LP solver)
    lambda_vec=logspace(-1,2,100); 
    for param=1:1:100
        lambda=lambda_vec(param);  
        sol=linprog(ones(2*m,1),[AA -AA; -AA AA],...
               [A'*b+ones(m,1)*sigma*lambda; -A'*b+ones(m,1)*sigma*lambda],...
               [],[],zeros(2*m,1),20*ones(2*m,1));   
        xDS=sol(1:m)-sol(m+1:2*m);
        Sout=find(abs(xDS)>1e-3);
        count=length(Sout)+1;
        disp(['    ==> lambda=',num2str(lambda),' and Count=',num2str(count)]);
        if count>50, break; end; 
        Er2_DSc{count,experiment}=[Er2_DSc{count,experiment},...
                                                   mean((xDS-x0).^2)/mean(x0.^2)]; 
        ErA2_DSc{count,experiment}=[ErA2_DSc{count,experiment},...
                                     mean((A*xDS-A*x0).^2)/mean((A*x0).^2)]; 
        xDS_proj=zeros(m,1); 
        xDS_proj(Sout)=pinv(A(:,Sout))*b;
        Er2_DS_projc{count,experiment}=[Er2_DS_projc{count,experiment},...
                                     mean((xDS_proj-x0).^2)/mean(x0.^2)]; 
        ErA2_DS_projc{count,experiment}=[ErA2_DS_projc{count,experiment},...
                                     mean((A*xDS_proj-A*x0).^2)/mean((A*x0).^2)]; 
    end;
    for count=1:1:50
        if ~isempty(Er2_DSc{count,experiment})
            Er2_DS(count,experiment)=min(Er2_DSc{count,experiment});
            ErA2_DS(count,experiment)=min(ErA2_DSc{count,experiment});
            Er2_DS_proj(count,experiment)=min(Er2_DS_projc{count,experiment});
            ErA2_DS_proj(count,experiment)=min(ErA2_DS_projc{count,experiment});
            Card_DS(count,experiment)=count-1;
        end;
    end;

    Best_BP=Best_BP+min(Er2_LARS_proj(find(Er2_LARS_proj(:,experiment)>0),experiment));
    Best_DS=Best_DS+min(Er2_DS_proj(find(Er2_DS_proj(:,experiment)>0),experiment));
    disp(['                           Best performance BP is ',num2str(Best_BP/experiment)]);
    disp(['                           Best performance DS is ',num2str(Best_DS/experiment)]);
    
    % Show the results
    figure(1); clf; 
    semilogy(0:1:49,mean(Er2_LARS(:,1:experiment),2),'b.'); 
    hold on;
    h=semilogy(0:1:49,mean(Er2_LARS_proj(:,1:experiment),2),'gs'); 
    set(h,'MarkerSize',4,'MarkerFaceColor','g');
    V=sum(Er2_DS(:,1:experiment),2)./sum(Er2_DS(:,1:experiment)>0,2);    
    h=semilogy(0:1:49,V,'ro'); 
    set(h,'MarkerSize',4,'MarkerFaceColor','r');
    V=sum(Er2_DS_proj(:,1:experiment),2)./sum(Er2_DS_proj(:,1:experiment)>0,2);    
    h=semilogy(0:1:49,V,'c^'); 
    set(h,'MarkerSize',6,'MarkerFaceColor','c');
    semilogy(0:50,ones(1,51)*mean(Er2_ORACLE(1:experiment))); 
    % axis([0 50 0.001 1]); 
    legend({'BP','BP-Projected','DS','DS-Projected','Oracle'},1);
    set(gca,'FontSize',14);
    h=xlabel('Cardinality: ||x||_0');
    set(h,'FontSize',14);
    h=ylabel('Accuracy: ||x-x_0||_2^2/||x_0||_2^2');
    set(h,'FontSize',14);

    figure(2); clf; 
    semilogy(0:1:49,mean(ErA2_LARS(:,1:experiment),2),'b.'); 
    hold on;
    h=semilogy(0:1:49,mean(ErA2_LARS_proj(:,1:experiment),2),'gs'); 
    set(h,'MarkerSize',4,'MarkerFaceColor','g');
    V=sum(ErA2_DS(:,1:experiment),2)./sum(ErA2_DS(:,1:experiment)>0,2);    
    h=semilogy(0:1:49,V,'ro'); 
    set(h,'MarkerSize',4,'MarkerFaceColor','r');
    V=sum(ErA2_DS_proj(:,1:experiment),2)./sum(ErA2_DS_proj(:,1:experiment)>0,2);    
    h=semilogy(0:1:49,V,'c^'); 
    set(h,'MarkerSize',6,'MarkerFaceColor','c');
    % axis([0 50 0.001 1]); 
    semilogy(0:1:50,ones(1,51)*mean(ErA2_ORACLE(1:experiment))); 
    legend({'BP','BP-Projected','DS','DS-Projected','Oracle'},1);
    set(gca,'FontSize',14);
    h=xlabel('Cardinality: ||x||_0');
    set(h,'FontSize',14);
    h=ylabel('Accuracy: ||Ax-Ax_0||_2^2/||Ax_0||_2^2');
    set(h,'FontSize',14);

    drawnow;    
   
end;

figure(1); clf;
h=semilogy(0:1:49,mean(Er2_LARS,2),'b');
set(h,'LineWidth',2);
hold on;
h=semilogy(0:1:49,mean(Er2_LARS_proj,2),'r');
set(h,'LineWidth',2);
%set(h,'MarkerSize',4,'MarkerFaceColor','g');
V=sum(Er2_DS,2)./sum(Er2_DS>0,2);
h=semilogy(0:1:49,V,'c');
set(h,'LineWidth',2);
% set(h,'MarkerSize',4,'MarkerFaceColor','r');
V=sum(Er2_DS_proj,2)./sum(Er2_DS_proj>0,2);
h=semilogy(0:1:49,V,'g');
set(h,'LineWidth',2);
% set(h,'MarkerSize',6,'MarkerFaceColor','c');
h=semilogy(0:50,ones(1,51)*mean(Er2_ORACLE(:,1)),'--');
set(h,'LineWidth',2);
axis([0 50 0.001 1]);
legend({'BP','BP-Projected','DS','DS-Projected','Oracle'},1);
set(gca,'FontSize',14);
h=xlabel('Cardinality: ||x||_0');
set(h,'FontSize',14);
h=ylabel('Accuracy: ||x-x_0||_2^2/||x_0||_2^2');
set(h,'FontSize',14);
% print -depsc2 Chapter_08_DS_vs_BP_rep.eps

figure(2); clf;
h=semilogy(0:1:49,mean(ErA2_LARS,2),'b');
set(h,'LineWidth',2);
hold on;
h=semilogy(0:1:49,mean(ErA2_LARS_proj,2),'r');
set(h,'LineWidth',2);
V=sum(ErA2_DS,2)./sum(ErA2_DS>0,2);
h=semilogy(0:1:49,V,'c');
set(h,'LineWidth',2);
V=sum(ErA2_DS_proj,2)./sum(ErA2_DS_proj>0,2);
h=semilogy(0:1:49,V,'g');
set(h,'LineWidth',2);
% axis([0 50 0.001 1]);
h=semilogy(0:1:50,ones(1,51)*mean(ErA2_ORACLE(:,1)),'--');
set(h,'LineWidth',2);
legend({'BP','BP-Projected','DS','DS-Projected','Oracle'},1);
set(gca,'FontSize',14);
h=xlabel('Cardinality: ||x||_0');
set(h,'FontSize',14);
h=ylabel('Accuracy: ||Ax-Ax_0||_2^2/||Ax_0||_2^2');
set(h,'FontSize',14);
axis([0 50 0.001 1]); 
% print -depsc2 Chapter_08_DS_vs_BP_sig.eps

return;

