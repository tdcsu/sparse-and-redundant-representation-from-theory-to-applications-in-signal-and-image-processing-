% Figure - NONE
% =========================================
% This code checks that the properties of the effective matrix A 
% (c and norms of atoms) are as found by the iterative methods

m=2^10; n=2^7; 
posH=randperm(m); posH=sort(posH(1:n));
Pin=spdiags((1:9/(m-1):10)',[0],m,m);
Pout=speye(m,m); Pout=Pout(posH,:);

% This builds the Hadamard matrix
H=zeros(m,m);
for k=1:1:m,
    disp(k/m);
    temp=zeros(m,1);
    temp(k)=1; 
    H(:,k)=fht(temp);
end;

A=Pout*H*Pin; 

% computing c for SSF using the iterative way
iter=100; 
temp=randn(m,1); 
for k=1:1:iter, 
    temp=temp/norm(temp); 
    temp=Pout*fht(Pin*temp);
    temp=Pin*ifht(Pout'*temp);
    disp(norm(temp));
end; 
c=norm(temp); 
disp(['The estimated c is ',num2str(c)]);
disp(['The TRUE c is ',num2str(norm(A'*A))]);


% computing diag(A'*A) for PCD using the iterative way
iter=1000; 
ww=zeros(m,1);
for k=1:1:iter,
    temp=randn(n,1); 
    temp=Pin*ifht(Pout'*temp);
    ww=ww+temp.^2;
end;
W=ww/iter;
figure(1); clf; 
plot(W); hold on; 
plot(diag(A'*A),'r'); 
legend({'Estimated norms','True norms'}); 

    