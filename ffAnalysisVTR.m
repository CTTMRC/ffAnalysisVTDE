function [vtde,euclideanDistance]=ffAnalysisVTR(X,Y,Delta,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expectedType=cellstr(["LS","GMM"]);
defaultOptimOptions=struct(...
                'nPop',300,...
                'maxGeneration',30,...
                'maxIter',30,...
                'CR',0.1,...
                'F',0.5...
                );   
P=inputParser;
addOptional(P,'VTRType',"LS",@(x)validCellString(x,expectedType));
addOptional(P,'Options',defaultOptimOptions,@(x)validOptionStructure(x,defaultOptimOptions));
addRequired(P,'Delta')
parse(P,Delta,varargin{:});
Delta=P.Results.Delta;
type=P.Results.VTRType;
optimOptions=P.Results.Options;


switch type
    case "LS"
        [vtde,euclideanDistance]=ffAnalysisVTRLS(X,Y,Delta,'Options',optimOptions);
    case "GMM"
        [vtde,euclideanDistance]=ffAnalysisVTRGMM(X,Y,Delta,'Options',optimOptions);
    otherwise
end
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
function OK=validOptionStructure(in,valid)
valid=fieldnames(valid);
if isstruct(in)
    idx=zeros(1,length(in));
    check=fields(in);
    for i=1:length(in)
        idx(i)=contains(check{i},valid);
    end
    OK=all(idx);
else
    OK=false;
end
end
