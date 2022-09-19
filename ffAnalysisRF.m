function [vtde,euclideanDistance]=ffAnalysisRF(X,Y,Delta)
P=inputParser;
addRequired(P,'Delta')
parse(P,Delta);
upper=P.Results.Delta.Upper;
lower=P.Results.Delta.Lower;
target=P.Results.Delta.Target;
maxLag=upper-lower;
span=[0,cumsum(maxLag)];
timeDelayIndex=zeros(size(maxLag));
maxMR=zeros(size(maxLag));
X=X{:,:};

treeTemplate = templateTree('Surrogate','on'...
    ,'MinLeafSize',30 ...
    ,'PredictorSelection','interaction-curvature');
regressionEnsemble = fitrensemble(X,Y...
    ,'Learners', treeTemplate...
    ,'Method','bag'...
    );

eOrig=loss(regressionEnsemble,X,Y);
eDivide=zeros(1,min(size(X)));
n_half=floor(max(size(X))/2);

eDivide(1)=loss(regressionEnsemble,[circshift(X(:,1),n_half),X(:,2:end)],Y);
for i=2:size(X,2)-1
    eDivide(i)=loss(regressionEnsemble,[X(:,1:i-1),circshift(X(:,i),n_half),X(:,i+1:end)]...
        ,Y);
end
eDivide(end)=loss(regressionEnsemble,[X(:,1:end-1),circshift(X(:,end),n_half)],Y);
estModelReliance=eDivide./eOrig;

for i=2:length(span)
    variableModelReliance=estModelReliance(span(i-1)+1:span(i));
    [maxMR(i-1),timeDelayIndex(i-1)]=max(abs(variableModelReliance));
end

% tic
% eOrig_0=zeros(1,min(size(X(:,:))));
% for i=1:max(size(X(:,:)))
%     eOrig_0(i)= loss(regressionEnsemble,X(i,:),Y(i,:));
% end
% eOrig_0=sum(eOrig_0)/max(size(X(:,:)));
%
% eDivide_0=zeros(1,min(size(X(:,:))));
% eDivide_temp=zeros(1,n_half);
% for j=1:n_half
%     eDivide_temp(j)= loss(regressionEnsemble,[X(n_half+j,1),X(j,2:end)],Y(j,:)) ...
%         + loss(regressionEnsemble,[X(j,1),X(n_half+j,2:end)],Y(n_half+j,:));
% end
% eDivide_0(1)=sum(eDivide_temp)/2/n_half;
% for i=2:269
%     eDivide_temp=zeros(1,n_half);
%     for j=1:n_half
%         eDivide_temp(j)= loss(regressionEnsemble,[X(j,1:i-1),X(n_half+j,i),X(j,i+1:end)],Y(j,:)) ...
%             + loss(regressionEnsemble,[X(n_half+j,1:i-1),X(j,i),X(n_half+j,i+1:end)],Y(n_half+j,:));
%     end
%     eDivide_0(i)=sum(eDivide_temp)/2/n_half;
% end
% eDivide_temp=zeros(1,n_half);
% for j=1:n_half
%     eDivide_temp(j)= loss(regressionEnsemble,[X(j,1:269),X(n_half+j,end)],Y(j,:))...
%         + loss(regressionEnsemble,[X(n_half+j,1:269),X(j,270)],Y(n_half+j,:));
% end
% eDivide_0(end)=sum(eDivide_temp)/2/n_half;
% ModelReliance=eDivide_0./eOrig_0;
% toc



vtde=max(timeDelayIndex-1,0);
distance=(vtde-(target-lower))./maxLag;
relevantDistance=distance([1:2,4:5,7:8]);
euclideanDistance=norm([relevantDistance,zeros(size(relevantDistance))]);
end
