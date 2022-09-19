function [vtde,euclideanDistance]=ffAnalysisPLS(X,Y,Delta,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING                                                                 %
% need to re write it without PLS toolbox for ease of access!!            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expectedType=cellstr(["SR","REG"]);
P=inputParser;
addOptional(P,'PLSType',"SR",@(x)validCellString(x,expectedType));
addRequired(P,'Delta')
parse(P,Delta,varargin{:});
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
type=P.Results.PLSType;
timeDelayIndex=zeros(size(maxLag));
maxIdx=zeros(size(maxLag));
X=X{:,:};
%WARNING THIS USES PLS TOOLBOX
plsOpt.preprocessing={'autoscale','autoscale'};
plsOpt.display='none';
plsOpt.plots='none';
tempMdl=pls(X,Y,max(maxLag),plsOpt);
tempMdl=crossval(X,Y,tempMdl,{'con' (6)},max(maxLag),plsOpt);
[~,nComp]=max(tempMdl.q2y);
plsMdl=pls(X,Y,nComp,plsOpt);
plsMdl=crossval(X,Y,plsMdl,{'con' (6)},max(maxLag),plsOpt);
switch type
    case "SR"
        for i=2:length(span)
            SR=plsMdl.selratio(span(i-1)+1:span(i));
            [maxIdx(i-1),timeDelayIndex(i-1)]=max(abs(SR));
        end
    case "REG"
        for i=2:length(span)
            REG=plsMdl.reg(span(i-1)+1:span(i));
            [maxIdx(i-1),timeDelayIndex(i-1)]=max(abs(REG));
        end
    otherwise
        
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
    OK=all(contains(in,valid));
end
end
