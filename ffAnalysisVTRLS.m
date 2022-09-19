function [vtde,euclideanDistance]=ffAnalysisVTRLS(X,Y,Delta,varargin)
defaultOptimOptions=struct(...
    'nPop',300,...
    'maxGeneration',30,...
    'maxIter',30,...
    'CR',0.1,...
    'F',0.5...
    );
P=inputParser;
addRequired(P,'Delta')
addOptional(P,'Options',defaultOptimOptions,@(x)validOptionStructure(x,defaultOptimOptions));
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
optimOptions=P.Results.Options;
nPop=optimOptions.nPop;
maxGeneration=optimOptions.maxGeneration;
maxIter=optimOptions.maxIter;
crossoverRate=optimOptions.CR;
scalingFactor=optimOptions.F;
checkDebug=fieldnames(optimOptions);
isDebug=sum(contains(checkDebug,"Debug"))>0;
if isDebug
    verbose=optimOptions.debug;
else
    verbose=false;
end
w=mean(X);
tau=1;
X=X';
vars=size(X,1);
obs=size(X,2);
Y=Y(max(maxLag) + 1 :end,:)';
for i = 1 : maxIter
    [vtde,trace]=IDE_VTRLS(vars,max(maxLag),nPop,maxGeneration,X,Y,w,tau,crossoverRate,scalingFactor);
    X_s = X(:, max(maxLag) + 1  : end ) ;
    for j = 1 : vars
        X_s(j, :) = circshift(X(j, max(maxLag) + 1  : end ) ,  vtde(j) );
    end
    w = (Y * X_s') / (X_s * X_s');
    tau = (1/(obs - max(maxLag))) * sum((Y - w * X_s) .* (Y - w * X_s));
end
if verbose
figure
bar(vtde)

MLE = -1 * (trace(:, 2));
figure
plot(MLE, '.-') 
end

distance=(vtde-(target-lower))./maxLag;
relevantDistance=distance([1:2,4:5,7:8]);
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);

end

function OK=validOptionStructure(in,valid)
valid=fields(valid);
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