function [X_lagged,Y_lagged,XY_lagged]=TTlagify_v2(X,Y,Delta)
[obs,var]=size(X);
dLower=Delta.Lower;
dUpper=Delta.Upper;
if min(dLower)<0
    lagYFlag=1;
    lagYAmount=-min(Delta.Lower);
    dUpper=dUpper-min(Delta.Lower);
    dLower=dLower-min(Delta.Lower);
    
else
    lagYFlag=0;
    lagYAmount=0;
end



max_lag=max(Delta.max_lag,max(dUpper));



if obs> max_lag
    X_lagged=timetable('Size',[obs,var*max_lag],'VariableTypes',repmat("doublenan",[1,var*max_lag]),'RowTimes',X.Time);
    Y_lagged=timetable('Size',[obs,size(Y,2)],'VariableTypes',repmat("doublenan",[1,size(Y,2)]),'RowTimes',Y.Time);
    Y_lagged.Properties.VariableNames=Y.Properties.VariableNames;
else
    warning(("The minimum number of observations is max_lag" + newline...
        + "Not enough observations: Aborting."));
    return
end

names_array=cell(var*max_lag,1);
Duration=mode(diff(X.Time));
intrvl=[0;find(diff(X.Time)~=Duration);obs];
k=1;
for l=1:length(intrvl)-1
    if (intrvl(l+1)-intrvl(l))>=max_lag
        X_Temp=X(intrvl(l)+1:intrvl(l+1),:);
        Y_Temp=Y(intrvl(l)+1:intrvl(l+1),:);
        if lagYFlag
            Y_partial=lag(Y_Temp,lagYAmount);
        else
            Y_partial=Y_Temp;
        end
        Y_lagged(intrvl(l)+1:intrvl(l+1),:)=Y_partial;
        for i=1:max_lag
            
            X_partial=lag(X_Temp,i-1);
            
            for j=1:var
                if i>dLower(j)&&i<=dUpper(j)
                    idx=max_lag*(j-1)+i;
                    X_lagged(intrvl(l)+1:intrvl(l+1),idx)=X_partial(:,j);
                    idx_names=max_lag*(j-1)+i;
                    if k==1
                        if lagYAmount==0
                            if i>1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i) + ")";
                                
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            end
                        else
                            if i>lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i-lagYAmount-1) + ")";
                                
                            elseif i==lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t+" + num2str(lagYAmount-i+1) + ")";
                            end
                        end
                        
                    end
                else
                    idx_names=max_lag*(j-1)+i;
                    if k==1
                        if lagYAmount==0
                            if i>1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i) + ")";
                                
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            end
                        else
                            if i>lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t-" + num2str(i-lagYAmount-1) + ")";
                                
                            elseif i==lagYAmount+1
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t)";
                            else
                                names_array{idx_names}=X.Properties.VariableNames{j} + " (t+" + num2str(lagYAmount-i+1) + ")";
                            end
                        end
                        
                    end
                end
            end
            
            
        end
        
        k=0;
    end
    
end

names_array=cellstr(names_array);
X_lagged.Properties.VariableNames=names_array;
X_lagged=X_lagged(max_lag+1:end,:);
Y_lagged=Y_lagged(max_lag+1:end,:);
filter0=~(any(ismissing(X_lagged),1));
filter1=~(any(ismissing(X_lagged(:,filter0)),2)|any(ismissing(Y_lagged),2));
filter2=~(any(isinf(X_lagged{:,:}),2)|any(isinf(Y_lagged{:,:}),2));
filter3=filter1&filter2;
X_lagged=X_lagged(filter3,:);
Y_lagged=Y_lagged(filter3,:);
XY_lagged=[X_lagged,Y_lagged];

%
% missing_variables=sum(sum(ismissing(X_lagged)))
end
