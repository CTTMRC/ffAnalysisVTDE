function [xCalib,ttXCalib,yCalib,xTest,ttXTest,yTest,Delta,LagVect]=ffAnalysisTEST(i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i=
% 1=1,1
% 2=2,1
% 3=3,1
% 4=1,2
% 5=2,2
% 6=3,2
% 7=3,1
% 8=3,2
% 9=3,3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ExperimentData,LagVect,Delta]=ffAnalysisINIT;
rng(1,'twister')
% LagVect=[45,40,35,30,25,20,15,10,5];
ffFactors=fullfact([3,3]);
sizeField=cellstr(["hundred","thousand","tenThousand"]);
complexityField=cellstr(["C1","C2","C3"]);
experiment=i;
xCalib=     ExperimentData.(sizeField{ffFactors(experiment,1)}).calib.X;
ttXCalib=   ExperimentData.(sizeField{ffFactors(experiment,1)}).calib.TT;
yCalib=     ExperimentData.(sizeField{ffFactors(experiment,1)}).calib...
.(complexityField{ffFactors(experiment,2)});
xTest=      ExperimentData.(sizeField{ffFactors(experiment,1)}).test.X;
ttXTest=    ExperimentData.(sizeField{ffFactors(experiment,1)}).test.TT;
yTest=      ExperimentData.(sizeField{ffFactors(experiment,1)}).test...
.(complexityField{ffFactors(experiment,2)});
ttXCalib=ttXCalib(:,(std(ttXCalib{:,:},1))>0);
ttXTest=ttXTest(:,(std(ttXTest{:,:},1))>0);