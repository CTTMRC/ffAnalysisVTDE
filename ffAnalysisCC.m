function [vtde,euclideanDistance]=ffAnalysisCC(X,Y,Delta,varargin)

expectedType=cellstr(["Pearson","Kendall","Spearman"]);
P=inputParser;
addOptional(P,'CCType',"Pearson",@(x)validCellString(x,expectedType));
addRequired(P,'Delta')
parse(P,Delta,varargin{:});
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
type=P.Results.CCType;
timeDelayIndex=zeros(size(maxLag));
maxCc=zeros(size(maxLag));
for i=2:length(span)
cc=corr(X{:,span(i-1)+1:span(i)},Y,"Type",type);
[maxCc(i-1),timeDelayIndex(i-1)]=max(abs(cc));
end
vtde=(timeDelayIndex-1);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=abs(distance([1:2,4:5,7:8]));
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end

function OK=validCellString(in,valid)

if iscellstr(in)||iscell(in)||isstring(in)
    idx=zeros(1,length(in));
    for i=1:length(in)
        idx(i)=contains(in{i},valid);
    end
    OK=all(idx);
else
    in=cellstr(in);
    OK=all(contain(in,valid));
end
end
