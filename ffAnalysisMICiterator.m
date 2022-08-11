function [vtde,euclideanDistance]=ffAnalysisMICiterator(X,Y,Delta,iter)

P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
timeDelayIndex=zeros(size(maxLag));
maxMic=zeros(size(maxLag));
inspectMic=[];
for i=2:length(span)
    xTemp=X{:,span(i-1)+1:span(i)};
    mic=zeros(1,maxLag(i-1));
    for j=1:size(xTemp,2)
        temp=mine(xTemp(:,j)',Y',0.6,iter);
        mic(j)=temp.mic;
    end
    inspectMic=[inspectMic,mic];
[maxMic(i-1),timeDelayIndex(i-1)]=max(abs(mic));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=pdist([relevantDistance;zeros(size(relevantDistance))],'euclidean');
% if iter==3
%     keyboard
% end
end