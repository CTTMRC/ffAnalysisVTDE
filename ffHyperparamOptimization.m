function bestHyperparameters=ffHyperparamOptimization(X,Y,max_lag,varargin)
%%ARGCHECK
defaultTreeN=300;
validScalarPosNum = @(x) isnumeric(x) && all(x > 0);
expectedHPO=cellstr(["all","MinLeafSize","NumPredictorstoSample","MaxNumSplits"]);
expectedOptType=cellstr(["directMAE","invertedMAE"]);
p=inputParser;
addRequired(p,'max_lag');
addOptional(p,'numTrees',defaultTreeN,validScalarPosNum);
addOptional(p,'HParamOpt',"all",@(x)validCellString(x,expectedHPO));
addOptional(p,'OptType','directMAE',@(x) any(validatestring(x,expectedOptType)));
parse(p,max_lag,varargin{:});
max_lag=p.Results.max_lag;
TreeN=p.Results.numTrees;
HParamOptRqst=cellstr(p.Results.HParamOpt);
if contains(HParamOptRqst,'all')
    HParamOptRqst=cellstr(["MinLeafSize","NumPredictorstoSample","MaxNumSplits"]);
end
OptType=validatestring(p.Results.OptType,expectedOptType);
Hyperparameters = [];
if contains(HParamOptRqst,'all')
    Hyperparameters = [optimizableVariable('MinLeafSize',[1,max_lag], 'Type','integer')...
        ;optimizableVariable('NumPredictorstoSample', [1,width(X)-1], 'Type','integer')...
        ;optimizableVariable('MaxNumSplits', [1,round(length(X)/max_lag)], 'Type','integer')...
        ];
else
    if contains(HParamOptRqst,'MinLeafSize')
        Hyperparameters = optimizableVariable('MinLeafSize',[1,max_lag], 'Type','integer');
    end
    if any(contains(HParamOptRqst,'NumPredictorstoSample'))&&isempty(Hyperparameters)
        Hyperparameters = optimizableVariable('NumPredictorstoSample',[1,max_lag], 'Type','integer');
    elseif any(contains(HParamOptRqst,'NumPredictorstoSample'))
        Hyperparameters = [Hyperparameters ...
            ;optimizableVariable('NumPredictorstoSample', [1,width(X)-1], 'Type','integer')];
    end
    if any(contains(HParamOptRqst,'MaxNumSplits'))&&isempty(Hyperparameters)
        Hyperparameters = optimizableVariable('MaxNumSplits',[1,round(length(X)/max_lag)], 'Type','integer');
    elseif any(contains(HParamOptRqst,'MaxNumSplits'))
        Hyperparameters = [Hyperparameters ...
            ;optimizableVariable('MaxNumSplits', [1,round(length(X)/max_lag)], 'Type','integer')];
    end
end
switch OptType
    case "invertedMAE"
        fun = @(hparams)Inverse_oobMAE(hparams,X,Y,TreeN);
    case "directMAE"
        fun = @(hparams)oobMAE(hparams,X,Y,TreeN);
    otherwise
        error('?????')
end
results = bayesopt(fun, Hyperparameters);
bestHyperparameters = bestPoint(results);
end

function MAE = Inverse_oobMAE(params,x,y,num_trees)
p=inputParser;
validScalarPosNum = @(x) isnumeric(x) && all(x > 0);
addRequired(p,'numTrees',validScalarPosNum)
parse(p,num_trees);
numTrees=p.Results.numTrees;
%%%ORDER NOT TO BE CHANGED USES BINARY NUMBER TO SELECT THE FUNCTION!!!
f={@forest3,@forest2,@forest23, @forest1,@forest13,@forest12,@forest123};
expectedHPO=cellstr(["MinLeafSize","NumPredictorstoSample","MaxNumSplits"]);
optTarget=params.Properties.VariableNames;
target=contains(expectedHPO,optTarget);
binaryTarget=num2str(target);
k=bin2dec(binaryTarget); %%THIS BINARY CODE
randomForest=f{k}(params,x,y,numTrees);
MAE = -1/oobQuantileError(randomForest);   % Uses the median by default.
end

function MAE = oobMAE(params,x,y,num_trees)
p=inputParser;
validScalarPosNum = @(x) isnumeric(x) && all(x > 0);
addRequired(p,'numTrees',validScalarPosNum)
parse(p,num_trees);
numTrees=p.Results.numTrees;
%%%ORDER NOT TO BE CHANGED USES BINARY NUMBER TO SELECT THE FUNCTION!!!
f={@forest3,@forest2,@forest23, @forest1,@forest13,@forest12,@forest123};
expectedHPO=cellstr(["MinLeafSize","NumPredictorstoSample","MaxNumSplits"]);
optTarget=params.Properties.VariableNames;
target=contains(expectedHPO,optTarget);
binaryTarget=num2str(target);
k=bin2dec(binaryTarget); %%THIS BINARY CODE
randomForest=f{k}(params,x,y,numTrees);
MAE = oobQuantileError(randomForest);     



% randomForest = "TreeBagger(numTrees,x,y,'Method','regression','OOBPrediction','on'";
% for i=1:size(optTarget,2)
% randomForest = randomForest +",'"+ optTarget{i} +"',params."+optTarget{i};
% end
% randomForest=randomForest +");";
% randomForest=eval(randomForest);
% MAE = oobQuantileError(randomForest);   
end
function OK=validCellString(in,valid)

if iscellstr(in)||iscell(in)
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
function randomForest=forest1(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'MinLeafSize', params.MinLeafSize...
    );
end
function randomForest=forest2(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'NumPredictorstoSample', params.NumPredictorstoSample...
    );
end
function randomForest=forest3(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'MaxNumSplits',params.MaxNumSplits...
    );
end
function randomForest=forest12(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'MinLeafSize', params.MinLeafSize...
    ,'NumPredictorstoSample', params.NumPredictorstoSample...
    );
end
function randomForest=forest13(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'MinLeafSize', params.MinLeafSize...
    ,'MaxNumSplits',params.MaxNumSplits...
    );
end
function randomForest=forest23(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'NumPredictorstoSample', params.NumPredictorstoSample...
    ,'MaxNumSplits',params.MaxNumSplits...
    );
end
function randomForest=forest123(params,x,y,numTrees)
randomForest = TreeBagger(numTrees,x,y,'Method','regression',...
    'OOBPrediction','on'...
    ,'MinLeafSize', params.MinLeafSize...
    ,'NumPredictorstoSample', params.NumPredictorstoSample...
    ,'MaxNumSplits',params.MaxNumSplits...
    );
end
