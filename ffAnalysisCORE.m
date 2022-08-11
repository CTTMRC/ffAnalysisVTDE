function [ResultsTruthfulness,ResultsR2]=ffAnalysisCORE()
[ExperimentData,LagVect,Delta]=ffAnalysisINIT;
rng(1)
% LagVect=[45,40,35,30,25,20,15,10,5];
ffFactors=fullfact([3,3]);
sizeFactors=cellstr(["100___","1000__","10000_"]);
complexityFactors=cellstr(["linear","interactions","power_ratio"]);
rowNames=cellstr(cat(2,cell2str(sizeFactors(ffFactors(:,1))'),cell2str(complexityFactors(ffFactors(:,2))')));
methods=cellstr([...
    "P_Correlation",...
    "S_Correlation",...
    "K_Correlation",...
    "MIC",...
    "N-norm",...
    "Procrustes",...
    "MI-Sousa",...
    "MI-Index",...
    "PLS-SR",...
    "PLS-T2",...
    "DistanceCorrelation",...
    "VTR-LS",...
    "VTR-COSE",...
    "GA-ITSS",...
    "RandomForest"]);
ResultsR2=table('Size',[9,15],'VariableTypes',repmat("doublenan",[1,15]),'VariableNames',methods,'RowNames',rowNames);
ResultsTruthfulness=table('Size',[9,15],'VariableTypes',repmat("doublenan",[1,15]),'VariableNames',methods,'RowNames',rowNames);
sizeField=cellstr(["fiveT","fiveH","fiveTH"]);
complexityField=cellstr(["C1","C2","C3"]);
tic
vtdePc=zeros(9,9);
vtdeSc=zeros(9,9);
vtdeKc=zeros(9,9);
vtdeMic=zeros(9,9);
vtdeMgx=zeros(9,9);
vtdeProcrustes=zeros(9,9);
vtdeMiSousa=zeros(9,9);
vtdeMiIndex=zeros(9,9);
vtdeDCorrelation=zeros(9,9);

for experiment=1:9
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
    [vtdePc(experiment,:),distancePC]=ffAnalysisCC(ttXCalib,yCalib,Delta,...
        "CCType","Pearson");
    ResultsTruthfulness{experiment,"P_Correlation"}=distancePC;
    [vtdeSc(experiment,:),distanceSC]=ffAnalysisCC(ttXCalib,yCalib,Delta,...
        "CCType","Spearman");
    ResultsTruthfulness{experiment,"S_Correlation"}=distanceSC;
    [vtdeKc(experiment,:),distanceKC]=ffAnalysisCC(ttXCalib,yCalib,Delta,...
        "CCType","Kendall");
    ResultsTruthfulness{experiment,"K_Correlation"}=distanceKC;
    [vtdeMic(experiment,:),distanceMIC]=ffAnalysisMIC(ttXCalib,yCalib,Delta);
    ResultsTruthfulness{experiment,"MIC"}=distanceMIC;
    [vtdeMgx(experiment,:),distanceMGX]=ffAnalysisMGX(xCalib,yCalib,Delta);
    ResultsTruthfulness{experiment,"N-norm"}=distanceMGX;
    [vtdeProcrustes(experiment,:),distancePRO]=ffAnalysisPROCRUSTES(ttXCalib,yCalib,Delta);
    ResultsTruthfulness{experiment,"Procrustes"}=distancePRO;
    [vtdeMiIndex(experiment,:),distanceMIIDX]=ffAnalysisMIINDEX(ttXCalib,yCalib,Delta);
    ResultsTruthfulness{experiment,"MI-Index"}=distanceMIIDX;
    [vtdeDCorrelation(experiment,:),distanceDCORR]=ffAnalysisDISTANCE(ttXCalib,yCalib,Delta);
    ResultsTruthfulness{experiment,"DistanceCorrelation"}=distanceDCORR;
    
    
    
end
toc