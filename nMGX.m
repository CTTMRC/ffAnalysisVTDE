function [N_j,d,mic]=nMGX(x,y)
max_lag=size(y,2);
N_j=zeros(1,max_lag);
mic=struct(...
    'mic',[]...
    ,'mas',[]...
    ,'mev',[]...
    ,'mcn' ,[]...
    ,'mcn_general' ,[]...
    ,'tic' ,[]...
    );

for i=1:max_lag
    [N_j(i),mic(i)]=MGX(x,y(:,i));
end
[~,d]=max(N_j);
end




function [N,z]=MGX(x,y)
n=zeros(1,4);
n(1)=corr(y,x,'Type','Pearson');
n(2)=corr(y,x,'Type','Kendall');
n(3)=corr(y,x,'Type','Spearman');
z=mine(y',x',0.6,15);
n(4)=z.mic;
n=abs(n);
N=norm(n.*n',2);
% N=norm(n,2);
end