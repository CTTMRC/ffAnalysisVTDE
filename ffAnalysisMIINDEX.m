function [vtde,euclideanDistance]=ffAnalysisMIINDEX(X,Y,Delta)
P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
timeDelayIndex=zeros(size(maxLag));
maxMi=zeros(size(maxLag));
inspectMi=[];
for i=2:length(span)
    xTemp=X{:,span(i-1)+1:span(i)};
    mi=zeros(1,size(xTemp,2));
    for j=1:size(xTemp,2)
        mi(j)=MIxnyn(xTemp(:,j),Y,2);
    end
%     inspectMi=[inspectMi,mi];
    [maxMi(i-1),timeDelayIndex(i-1)]=max(abs(mi));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
