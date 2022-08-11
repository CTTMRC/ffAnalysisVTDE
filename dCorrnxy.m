function dCorr=dCorrnxy(x,y)

SF_Y=squareform(pdist(y,'euclidean'));
Bbar=repmat((mean(mean(SF_Y))),size(SF_Y));
Bk=repmat(mean(SF_Y,2),[1,size(SF_Y,2)]);
Bl=repmat(mean(SF_Y,1),[size(SF_Y,1),1]);
Bkl=SF_Y -Bl -Bk +Bbar;
V2Y=mean(mean(Bkl.*Bkl));
SF_X=squareform(pdist(x,'euclidean'));
Abar=repmat((mean(mean(SF_X))),size(SF_X));
Ak=repmat(mean(SF_X,2),[1,size(SF_X,2)]);
Al=repmat(mean(SF_X,1),[size(SF_X,1),1]);
Akl=SF_X -Al -Ak +Abar;
V2XY=mean(mean(Akl.*Bkl));
V2X=mean(mean(Akl.*Akl));
dCorr=sqrt(abs(V2XY/sqrt(V2X*V2Y)));
end
