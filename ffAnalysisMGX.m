function [vtde,euclideanDistance]=ffAnalysisMGX(X,Y,Delta)

P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
target=target-min(lower);
upper=upper-min(lower);
lower=lower-min(lower);
maxLag=upper-lower;
timeDelayIndex=zeros(size(maxLag));
yLag=max(upper);
yAugmented=repmat(Y,[1,yLag]);
NormJ=zeros(size(X,2),yLag);
for j=1:(yLag-1)
   yAugmented(:,j+1)=circshift(Y,-j);
end
for i=1:size(X,2)
    [NormJ(i,lower(i)+1:upper(i)),timeDelayIndex(i)]=nMGX(X(:,i),yAugmented(:,lower(i)+1:upper(i)));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
