function [vtde,euclideanDistance]=ffAnalysisGAITSS(X,Y,Delta)
P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
timeDelayIndex=zeros(size(maxLag));
s0=zeros(size(X));
sMax=0;
flag=0;
matInd=1:size(X,2);
inclusion=0;
iter=0;
includedVars=1:size(X,2);
X=X{:,:};
while flag==0
    iter=iter+1;
    MIxy=zeros(1,size(X,2));
    
    for i=1:size(X,2)
  
        MI=MIxnyn([s0(:,sum(s0,1)>0),X(:,i)],Y,6);%- MIxnyn(s0(:,sum(s0,1)>0),X(:,i),6);
       
        MIxy(:,i)=MI;
        
         
    end

    [~,j]=max((MIxy));
    
    if MIxy(j)>sMax
        
        includedVars=setdiff(includedVars,includedVars(j));
        inclusion=inclusion+1;
        sMax=MIxy(j);
        s0(:,inclusion)=X(:,j);
        X=X(:,setdiff(1:size(X,2),j));
        
        inclVAR=find((matInd(j)./(span+1))<=1,1,'first' )-1;
        inclVAR=max(inclVAR,1);
        inclLAG=mod(matInd(j),span(inclVAR));
        
        if inclLAG==0
            inclLAG=1;
        end
        
        matInd=matInd(setdiff(1:size(X,2)+1,j));
        timeDelayIndex(inclVAR)=inclLAG;
    else

        flag=1;
    end
    
end


vtde=max(timeDelayIndex-1,0);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=distance([1:2,4:5,7:8]);
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
