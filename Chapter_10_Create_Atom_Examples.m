% Figure - 10.4
% =========================================
% This program generates portions of Figure 10.4 - generating random 
% atoms from the redundant Haar dictionary.


function []=Chapter_10_Create_Atom_Examples(Number)

% This function generates a the redundant Haar matrix of size 400*2800, 
% representing 3 un-decimated levels of resolution. Then it generates a
% set of "Number" (anything between 1 and 99 should be fine) atoms, 
% extracted from it at random

d=20;
D1=zeros(d);
v=[1 zeros(1,d-2), -1];
for k=1:1:d
    D1(k,:)=v;
    v=[v(end),v(1:end-1)];
end;
D2=zeros(d);
v=[1 1 zeros(1,d-4), -1 -1];
for k=1:1:d
    D2(k,:)=v;
    v=[v(end),v(1:end-1)];
end;

S1=abs(D1);
S2=abs(D2);
Haar=[kron(S2,S2),kron(S2,D2),kron(D2,S2),kron(D2,D2),...
                             kron(S1,D1),kron(D1,S1),kron(D1,D1)];
for k=1:1:7*d^2
    Haar(:,k)=Haar(:,k)/sum(abs(Haar(:,k)));
end;

for j=11:1:Number+10 
    k=round(rand*2800);
    V=Haar(:,k); 
    V=V/max(V); 
    V=V*100+127;
    V=reshape(V,[20,20]);
    image(V); axis image; axis off; 
    colormap(gray(256)); 
    % uncomment the lines below to get the examples as eps files
    % FileName = ['Atom_Example_',num2str(j),'.eps'];
    % eval(['print -depsc2 ',FileName]);
    pause(0.1); 
end;

return;

% The following code generates the complete transform matrix in an
% laternative way 
% d=20^2; % size of the image
% n=d*7; % number of atoms
% H=zeros(400,2800);
% for k=0:1:2799 
%     H(:,k+1)=WaveletFiguBlur0(1,d,n,[zeros(k,1); 1; zeros(n-k-1,1)],1:1:n, n); 
% end;
