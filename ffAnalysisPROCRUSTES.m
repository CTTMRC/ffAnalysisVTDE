function [vtde,euclideanDistance]=ffAnalysisPROCRUSTES(X,Y,Delta)
P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
max_lag=upper-lower;
span=[0,cumsum(max_lag)];
timeDelayIndex=zeros(size(max_lag));
ProcrustesMax=zeros(size(max_lag));
Y=normalize(Y);
for i=2:length(span)
    xTemp=X{:,span(i-1)+1:span(i)};
    xWhole=reshape(xTemp,1,[]);
    xTemp=(xTemp-mean(xWhole))./std(xWhole);
    ProcrustesDistance=zeros(size(xTemp,2),1);
    for j=1:size(xTemp,2)
        ProcrustesDistance(j)=procrustes(xTemp(:,j),Y);
    end
[ProcrustesMax(i-1),timeDelayIndex(i-1)]=min(abs(ProcrustesDistance));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./max_lag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
