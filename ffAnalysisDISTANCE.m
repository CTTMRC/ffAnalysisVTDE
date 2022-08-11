function [vtde,euclideanDistance]=ffAnalysisDISTANCE(X,Y,Delta)
P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
max_lag=upper-lower;
span=[0,cumsum(max_lag)];
timeDelayIndex=zeros(size(max_lag));
dCorrelationMax=zeros(size(max_lag));
Y=normalize(Y);
for i=2:length(span)
    xTemp=X{:,span(i-1)+1:span(i)};
    xWhole=reshape(xTemp,1,[]);
    xTemp=(xTemp-mean(xWhole))./std(xWhole);
    dCorrelationDistance=zeros(size(xTemp,2),1);
    for j=1:size(xTemp,2)
        dCorrelationDistance(j)=dCorrnxy(xTemp(:,j),Y);
    end
[dCorrelationMax(i-1),timeDelayIndex(i-1)]=min(abs(dCorrelationDistance));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./max_lag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
