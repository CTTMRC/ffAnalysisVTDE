%FULLFACTORIAL [size + complexity]
function [ExperimentData,LagVector,Delta]=ffAnalysisINIT

autoRegressionCoefficients=[0.9,0.9,0.1];
% dLower=[30,25,20,15,10,5,0,-5,-10];
% dUpper=[60,55,50,45,40,35,30,25,20];
dLower=[0,0,0,0,0,0,0,0,0];
dUpper=[30,30,30,30,30,30,30,30,30];
settingOutliers=true;
settingStartSet=autoRegressionCoefficients;
settingYNoise=5;
settingLagNumber=3;
settingXNoise=[0.1,0.1,0.1];
rng(1)
[xSize1Calib,~,ttXSize1Calib,ttYSize1Calib,LagVector,~,~,Delta]=SynthDataGeneration(100,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );
[xSize2Calib,~,ttXSize2Calib,ttYSize2Calib]=SynthDataGeneration(1000,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );
[xSize3Calib,~,ttXSize3Calib,ttYSize3Calib]=SynthDataGeneration(10000,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );
rng('shuffle')
[xSize1Test,~,ttXSize1Test,ttYSize1Test]=SynthDataGeneration(100,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );
[xSize2Test,~,ttXSize2Test,ttYSize2Test]=SynthDataGeneration(1000,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );
[xSize3Test,~,ttXSize3Test,ttYSize3Test]=SynthDataGeneration(10000,1,...
    'Outliers',settingOutliers,...
    'startSet',settingStartSet,...
    'YNoise',settingYNoise,...
    'lagNum',settingLagNumber,...
    'XNoise',settingXNoise,...
    'deltaUpper',dUpper,...
    'deltaLower',dLower...
    );

%%%
%calib
ExperimentData.fiveT.calib.X=xSize1Calib;
ExperimentData.fiveT.calib.TT=ttXSize1Calib;
ExperimentData.fiveT.calib.C1=ttYSize1Calib.y_lin;
ExperimentData.fiveT.calib.C2=ttYSize1Calib.y_interactions;
ExperimentData.fiveT.calib.C3=ttYSize1Calib.y_enzyme;
%test
ExperimentData.fiveT.test.X=xSize1Test;
ExperimentData.fiveT.test.TT=ttXSize1Test;
ExperimentData.fiveT.test.C1=ttYSize1Test.y_lin;
ExperimentData.fiveT.test.C2=ttYSize1Test.y_interactions;
ExperimentData.fiveT.test.C3=ttYSize1Test.y_enzyme;
%%%
%calib
ExperimentData.fiveH.calib.X=xSize2Calib;
ExperimentData.fiveH.calib.TT=ttXSize2Calib;
ExperimentData.fiveH.calib.C1=ttYSize2Calib.y_lin;
ExperimentData.fiveH.calib.C2=ttYSize2Calib.y_interactions;
ExperimentData.fiveH.calib.C3=ttYSize2Calib.y_enzyme;
%test
ExperimentData.fiveH.test.X=xSize2Test;
ExperimentData.fiveH.test.TT=ttXSize2Test;
ExperimentData.fiveH.test.C1=ttYSize2Test.y_lin;
ExperimentData.fiveH.test.C2=ttYSize2Test.y_interactions;
ExperimentData.fiveH.test.C3=ttYSize2Test.y_enzyme;
%%%
%calib
ExperimentData.fiveTH.calib.X=xSize3Calib;
ExperimentData.fiveTH.calib.TT=ttXSize3Calib;
ExperimentData.fiveTH.calib.C1=ttYSize3Calib.y_lin;
ExperimentData.fiveTH.calib.C2=ttYSize3Calib.y_interactions;
ExperimentData.fiveTH.calib.C3=ttYSize3Calib.y_enzyme;
%test
ExperimentData.fiveTH.test.X=xSize3Test;
ExperimentData.fiveTH.test.TT=ttXSize3Test;
ExperimentData.fiveTH.test.C1=ttYSize3Test.y_lin;
ExperimentData.fiveTH.test.C2=ttYSize3Test.y_interactions;
ExperimentData.fiveTH.test.C3=ttYSize3Test.y_enzyme;
%%%






