function [D,G,Res]=DesignGrass(N,L,Iter,dd1,dd2,Init)
%========================================
% In this program we iteratively build a frame that have 
% the following properties:
%   1. normalized columns
%   2. minimal M (mutual incoherence)
% 
% Inputs:    N - Number of rows in the dictionary D
%               L - Number of columns in the dictionary D
%               Iter - Number of iterations
%               dd1 - Relative number of the G entries to shrink
%               dd2 - Shrink factor to use (1/dd2 is used to expand)
%               Init - An initialization dictionary
%
% Outputs:  D - The dictionary
%               G - The Gram matrix
%
% Example: [D,G]=DesignFrame(100,500,100,0.8,0.9);
%========================================

if nargin==5,
    D=randn(N,L); % initialization
    D=D*diag(1./sqrt(diag(D'*D))); % normalize columns
else, 
    D=Init;
end;
G=D'*D; % compute the Gram matrix
mu=sqrt((L-N)/N/(L-1)); 

Res=zeros(Iter,3); 
for k=1:1:Iter,
    
    % shrink the high inner products
    gg=sort(abs(G(:))); 
    pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    
    % expand the near zero products
    % pos=find(abs(G(:))<gg(round((1-dd1)*(L^2-L))));
    % G(pos)=G(pos)/dd2;

    % reduce the rank back to N
    [U,S,V]=svd(G); 
    S(N+1:end,1+N:end)=0;
    G=U*S*V';
    
    % Normalize the columns
    G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G)));

    % Show status
    gg=sort(abs(G(:))); 
    pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    fprintf(1,'%6i  %12.8f   %12.8f  %12.8f \n',...
        [k,mu,mean(abs(G(pos))),max(abs(G(pos)))]);
    
    Res(k,:)=[mu,mean(abs(G(pos))),max(abs(G(pos)))];

end;

[U,S,V]=svd(G); 
D=S(1:N,1:N).^0.5*U(:,1:N)';

return;



