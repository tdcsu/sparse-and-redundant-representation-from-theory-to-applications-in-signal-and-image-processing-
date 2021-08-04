% Figure - NONE
% =========================================
% In this script we explore the efficiency of the weak-MP, 
% compared to the plain MP. The efficiency is given as a ratio in the 
% range [0,1], repesentating the relative number of inner-products
% performed by the weak-MP, compared to the MP. 

n=30; m=50; Smax=10; Exper=2000; 

A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 

Efficiency=zeros(Smax,1);

count=zeros(Smax,Exper); 
for S=1:1:Smax,
   
    for experiment=1:1:Exper
       
        % Generate a test signal of cardinality S
        x=zeros(m,1);
        pos=randperm(m);
        x(pos(1:S))=sign(randn(S,1)).*(1+rand(S,1));
        b=A*x; 
    
        % Apply WMP
        thrWMP=1e-8; t=0.5; 
        r=b;
        xWMP=zeros(m,1);
        for kk=1:1:S,
            Z=abs(A'*r);
            posZ=find(Z>=t*sqrt(r'*r),1);
            if isempty(posZ)
                count(S,experiment)=count(S,experiment)+1/S; 
                posZ=find(Z==max(Z),1);
            else
                count(S,experiment)=count(S,experiment)+posZ/m/S; 
            end;
            xWMP(posZ)=xWMP(posZ)+A(:,posZ)'*r;
            r=r-A(:,posZ)*A(:,posZ)'*r;
        end;
        
    end;    
    Efficiency(S)=mean(count(S,:)); 
end;

figure(1); clf; 
h=plot(1:1:Smax,Efficiency,'b');  
set(h,'LineWidth',2);
h=xlabel('Cardinality of the true solution'); set(h,'FontSize',14);
h=ylabel('Efficiency Factor'); set(h,'FontSize',14);
set(gca,'FontSize',14);
axis([1 Smax 0 1]); 





