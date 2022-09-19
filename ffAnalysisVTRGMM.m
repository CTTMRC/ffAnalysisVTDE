function [vtde,euclideanDistance]=ffAnalysisVTRGMM(X,Y,Delta,varargin)
defaultOptimOptions=struct(...
    'nPop',300,...
    'maxGeneration',50,...
    'maxIter',50,...
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

X=X';
Y=Y(max(maxLag) + 1 :end,:)';
xVars=size(X,1);
vars=size(X,1)+size(Y,1);
obs=size(X,2);

K = 1;

% priors
alpha0 = ones(K,1);
m0 = zeros(1,vars);
beta0 = 1;
W0 = eye(vars)*0.1;
v0 = max(maxLag);

% initialize variational hyper-parameters
Pi = rand(1,K);
Pi=Pi/sum(Pi+eps);
beta = ones(K,1);
W=repmat(eye(vars) * 100,[1,1,K]);
m = randn(K, vars);
v = repmat(v0,K,1);
L=zeros(1,maxIter);
tic
rng(1)
for iter=1:maxIter
        
    [vtde,~]=IDE_VTRGMM(xVars, max(maxLag), nPop, maxGeneration, X, Y, Pi, beta, m, v, W, alpha0, m0, beta0, W0, v0,crossoverRate,scalingFactor);
    X_S = X(:, max(maxLag) + 1  : end ) ;
    for i = 1 : xVars
        X_S(i, :) = X(i, max(maxLag) + 1 - vtde(i) : obs - vtde(i));
    end
    X0 = [X_S', Y'];
    phi=update_phi(Pi, m, v, W, beta, X0, K, vars);
    Nk=sum(phi, 1)';
    Pi=alpha0+Nk;
    beta=beta0+Nk;
    v=v0+Nk;
    m=update_m(phi, beta, m0, beta0, X0, vars);
    W=update_w(phi, beta, m, W0, beta0, m0, X0, vars);
    %ELBO
    lb = elbo(phi, Pi, beta, v, W, alpha0, beta0, v0, W0);
    L(iter) = lb;
   
end
if verbose
    figure
    bar(vtde)
    
    figure(2)
    plot(L,'b.-');
    Spred=zeros(size(W));
    for i = 1: K
        Spred(:, :, i) = W(:, :, i) / v(i, :);
    end
end


toc
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