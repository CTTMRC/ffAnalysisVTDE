function [x,y_lin,ttX,ttY,LagVector,Coeff,InnerX,Delta]=SynthDataGeneration(n,t,varargin)
%%ARGCHECK
defaultSet = [0;0;0];
defaultFull = repmat([0;0;0],1,t);
defaultDeltaLower=zeros(1,9);
defaultDeltaUpper=30*ones(1,9);
defaultPartitionNumber = 1;
defaultMultinomial = [0.75/6,0.75/6,0.75/6,0.25,0.75/6,0.75/6,0.75/6] ;
defaultXXXnoise=[1,1,1];
defaultYnoise=0.1;
defaultLag=3;
P = inputParser;
validScalarNonnegNum = @(x) isnumeric(x) && all(x >= 0);
validScalarPosNum = @(x) isnumeric(x) && all(x > 0);
validFullSet= @(x) validateCoeffcients(x) && all(size(x)==[3,t]);
validMultinomial= @(x) isnumeric(x) && all(x >= 0) && sum(x)==1 && mod(numel(x),2)==1;
addOptional(P,'fullSet',defaultFull,validFullSet)
addOptional(P,'deltaUpper',defaultDeltaUpper,@(x) isnumeric(x))
addOptional(P,'deltaLower',defaultDeltaLower,@(x) isnumeric(x))
addOptional(P,'Xnoise',defaultXXXnoise,validScalarNonnegNum)
addOptional(P,'Ynoise',defaultYnoise,validScalarNonnegNum)
addOptional(P,'Outliers',false)
addOptional(P,'startSet',defaultSet,validScalarNonnegNum)
addOptional(P,'Verbose',false)
addOptional(P,'isDynamic',false)
addOptional(P,'DynamicPartitions', defaultPartitionNumber, validScalarPosNum)
addOptional(P,'lagNum', defaultLag, validScalarPosNum)
addOptional(P,'MultinomialProbabilities', defaultMultinomial, validMultinomial)
parse(P,varargin{:});
c0=P.Results.startSet;
if length(c0)==1
    c0=[c0;c0;c0];
end
fullSet=P.Results.fullSet;
Verbose=P.Results.Verbose;
Outliers=P.Results.Outliers;
isDynamic=P.Results.isDynamic;
DynamicPartitions=P.Results.DynamicPartitions;
MNProb=P.Results.MultinomialProbabilities;
NoiseY=P.Results.Ynoise;
NoiseX=P.Results.Xnoise;
Delta.Lower=P.Results.deltaLower;
Delta.Upper=P.Results.deltaUpper;
lag=P.Results.lagNum;
%%INIT_PARAM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expressedVariables=9;
LagVector=repmat(fliplr(lag*(1:expressedVariables)),[DynamicPartitions,1]);
lag_max_lead=lag+lag*expressedVariables;
Delta.max_lag=lag_max_lead;
Delta.Target=LagVector;
t_safety=101*t;
safety_tail=max(101*t,lag_max_lead+1);
k=n+t_safety+safety_tail;
coeff_orig=[4.123,-4.234,4.345];
coeff_orig=coeff_orig./norm(coeff_orig);
x1_coeff=circshift(coeff_orig,0); %[ 0.5621  -0.5772   0.5924];
x2_coeff=circshift(coeff_orig,1); %[ 0.5924   0.5621  -0.5772];
x3_coeff=circshift(coeff_orig,2); %[-0.5772   0.5924   0.5621];
lin_coeff=[-0.7071,0.7071]; %Norm==1
int_coeff=[-0.2300    0.2300   -0.7146   -0.4732    0.3995];%Norm==1
safe_n=n+safety_tail;
kernel_coeff=((1:safe_n).^2.*sin((1:safe_n).^2))./(1+(1:safe_n).*cos((1:safe_n)));
kernel_coeff=kernel_coeff/norm(kernel_coeff);
sigma_kernel=37/23;
gamma_kernel=1/(2*(sigma_kernel)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%MODELS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linear=@(X)normalize(X(:,1)*lin_coeff(1)...
    +X(:,2)*lin_coeff(2),'medianiqr');
interactions=@(X)normalize(X(:,1)*int_coeff(1)...
    +X(:,2)*int_coeff(2)...
    +X(:,1).^2*int_coeff(3)...
    +X(:,2).^2*int_coeff(4)...
    +X(:,1).*X(:,2)*int_coeff(5),'medianiqr');
nonLinear=@(X)normalize(sign((X(:,1))*lin_coeff(1)).*log2(abs(X(:,1)*lin_coeff(1)))...
    +sign((X(:,2))*lin_coeff(2)).*nthroot(abs(X(:,2)*lin_coeff(2)),3),'medianiqr');
kernel=@(X)normalize(X*kernel_coeff','medianiqr');
exponential=@(X)normalize(exp(nthroot(X(:,1)*lin_coeff(1),5).^2 ...
    +nthroot(X(:,2)*lin_coeff(2),5).^2));
powerRatio=@(X)normalize(power(((X(:,2))*lin_coeff(2)),2)...
./(1+power(((X(:,1)).*lin_coeff(1)),3)),'medianiqr');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filter_stable = false;
c1=zeros(t,1);
c2=zeros(t,1);
c3=zeros(t,1);
if all(fullSet==defaultSet)
    count=0;
    while ~filter_stable
        count=count+1;
        if all(c0==defaultSet)
            c0=rand(3,1);
        end
        c1(1)=c0(1);
        c2(1)=c0(2);
        c3(1)=c0(3);
        c1(2:t)=randperm(100,t-1)*0.001*(1-c1(1));%power(c1(1),[1,10*randperm(10*t,t-1)]);
        c2(2:t)=randperm(100,t-1)*0.001*(1-c1(1));%power(c2(1),[1,10*randperm(10*t,t-1)]);
        c3(2:t)=randperm(100,t-1)*0.001*(1-c1(1));%power(c3(1),[1,10*randperm(10*t,t-1)]);
        coefficients1 = [1; -c1];
        coefficients2 = [1; -c2];
        coefficients3 = [1; -c3];
        if max(abs(roots(coefficients1))) < 1 && max(abs(roots(coefficients2))) < 1 && max(abs(roots(coefficients3))) < 1
            filter_stable = true;
            b1=c1;
            b2=c2;
            b3=c3;
        end
        if count >n*t*1000
            msg="No stable filter configurations."...
                +newline+"Try:"...
                +newline+" - Different starting points"...
                +newline+" - Different seed"...
                +newline+" - Lower autoregression order"...
                ;
            error(msg)
        end
    end
else
    b1=fullSet(1,:);
    b2=fullSet(2,:);
    b3=fullSet(3,:);
end
Coeff=[b1';b2';b3'];

if Verbose>=1
    CoeffStr="";
    for i=1:3
        CoeffStr=CoeffStr+string(i)+": ";
        for j=1:t
            CoeffStr=CoeffStr+string(Coeff(i,j))+", ";
        end
        CoeffStr=CoeffStr+newline;
    end
    msg="Generating data with AR order "+string(t)...
        + newline +CoeffStr;
    disp(msg);
end
base=mvnrnd(zeros(3,1),eye(3),k);
X1=filter(1,[1;-b1],base(:,1));
X2=filter(1,[1;-b2],base(:,2));
X3=filter(1,[1;-b3],base(:,3));
X1=X1(t_safety+1:end,:);
X2=X2(t_safety+1:end,:);
X3=X3(t_safety+1:end,:);
X1=normalize(X1,'zscore');
X2=normalize(X2,'zscore');
X3=normalize(X3,'zscore');
x=zeros(size(X2,1),9);
X1_x=X1*x1_coeff;
X2_x=X2*x2_coeff;
X3_x=X3*x3_coeff;
if Verbose
    disp(char(hex2dec('2713')))
end
% creation of matrix X
for i=1:3
    XNoise=mvnrnd(zeros(3,1),diag(NoiseX),size(x,1));
    x1=(i*3-2);
    x2=(i*3-1);
    x3=(i*3);
    %"Measurement" error on the measured x
    x(:,x1)=X1_x(:,i)+XNoise(:,1);
    x(:,x2)=X2_x(:,i)+XNoise(:,2);
    x(:,x3)=X3_x(:,i)+XNoise(:,3);
end
%%INSERT OUTLIERS

if Outliers
    if Verbose
        msg="Adding Outliers...";
        disp(msg)
        
    end
    OUTProb1=0.05;
    OUTProb2=[3/6,2/6,1/6];
    OuliersCount=0;
    span_outliers=lag_max_lead+1:size(x,1)-(safety_tail-lag_max_lead);
    for i=span_outliers
        r=binornd(1,OUTProb1);
        if r
            OuliersCount=OuliersCount+1;
            r2=(1:3).*mnrnd(1,OUTProb2,1);
            r2(r2==0)=[];
            if  r2==3
                P=randperm(9,3);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*1.*iqr(x(:,P));
                
            elseif r2==2
                P=randperm(9,2);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*2.*iqr(x(:,P));
                
            elseif r2==1
                P=randperm(9,1);
                x(i,P)=x(i,P)+sign(x(i,P)-median(x(:,P)))*3*iqr(x(:,P));
                
            end
        end
    end
    totalOutliers=(OuliersCount/n)*100;
    if Verbose
        msg=string(sprintf('%2.1f',round(totalOutliers,1)))+"% of total observations";
        disp(msg)
        disp(char(hex2dec('2713')))
    end
end

%%INSERT LAG AND DYNAMIC SRUCTURE
switch isDynamic
    case false
        j=1;
        for i=fliplr(1:9)
            VariableLag=-lag*j;
            x(:,i)=circshift(x(:,i),VariableLag);
            j=j+1;
        end
    case true
        Dyn_statesNum=numel(MNProb);
        Dyn_range=floor(Dyn_statesNum/2);
        Dyn_states=-Dyn_range:Dyn_range;
        Dyn_size=round(n/DynamicPartitions);
        span_Dynamic=lag_max_lead+1:size(x,1)-(safety_tail-lag_max_lead);
        l=1;
        if Verbose
            msg="Adding Dynamic Delay..."+newline+"Delay Changes every "+string(Dyn_size)+" points.";
            disp(msg)
        end
        
        for i=fliplr(1:9)
            VariableLag=-lag*l;
            x(:,i)=circshift(x(:,i),VariableLag);
            LaggedVariable=x(span_Dynamic,i);
            if DynamicPartitions >1
                for j=2:DynamicPartitions
                    lower_x=span_Dynamic(1)+((j-1)*Dyn_size);
                    upper_x=span_Dynamic(1)+((j)*Dyn_size)-1;
                    lower_lag=(j-1)*Dyn_size+1;
                    upper_lag=j*Dyn_size;
                    r=mnrnd(1,MNProb,1);
                    R=Dyn_states(r==1);
                    Temp_Lag=circshift(LaggedVariable,R);
                    x(lower_x:upper_x,i)=Temp_Lag(lower_lag:upper_lag);
                    LagVector(j,i)=LagVector(j,i)-R;
                end
            else
            end
            l=l+1;
        end
        if Verbose
            
            disp(char(hex2dec('2713')))
        end
end
%%CREATE Y
D=([X1,X2]);%normalize
yNoise=NoiseY*mvnrnd(zeros(6,1),eye(6),size(D,1));
DD=squareform(pdist(D,'squaredeuclidean'));
Gauss_kernel=exp(-gamma_kernel*sqrt(DD));
y_lin=linear(D)+yNoise(:,1);
y_int=interactions(D)+yNoise(:,2);
y_nlin=nonLinear(D)+yNoise(:,3);
y_kernel=kernel(Gauss_kernel)+yNoise(:,4);
y_exp=exponential(D)+yNoise(:,5);
y_enzyme= powerRatio(D)+yNoise(:,6);
%%CREATE TABLES
time_idx=datetime('yesterday')+minutes(0:lag:lag*(size(D,1)-1));
ttX = timetable('Size',[size(D,1),9],'VariableTypes',repmat("double",[1,9]),'RowTimes',time_idx);
ttY = timetable('Size',[size(D,1),6],'VariableTypes',repmat("double",[1,6]),'RowTimes',time_idx);
ttX.Properties.VariableNames={'X1_1','X2_1','X3_1','X1_2','X2_2','X3_2','X1_3','X2_3','X3_3'};
ttY.Properties.VariableNames={'y_lin','y_interactions','y_kernel','y_nlin','y_exp','y_enzyme'};
ttX{:,:}=x;
ttY{:,:}=[y_lin,y_int,y_kernel,y_nlin,y_exp,y_enzyme];
[ttX,ttY]=TTlagify_v2(ttX,ttY,Delta);
if  max(Delta.Upper)>lag_max_lead
    lag_max_lead=max(Delta.Upper);
end
if  min(Delta.Lower)<0
    lag_max_lead=lag_max_lead-min(Delta.Lower);
end
Delta.max_lag=lag_max_lead;
% [ttX,ttY]=TTlagify(ttX,ttY,lag_max_lead);
x=x(lag_max_lead+1:end-(safety_tail-lag_max_lead),:);
y_lin=y_lin(lag_max_lead+1:end-(safety_tail-lag_max_lead),:);
ttY=ttY(1:end-(safety_tail-lag_max_lead),:);
ttX=ttX(1:end-(safety_tail-lag_max_lead),:);

InnerX=[X1,X2,X3];
InnerX=InnerX(lag_max_lead+1:end-(safety_tail-lag_max_lead),:);
end


% %%OTHER_KERNELS
% Gauss_kernel=exp(-sqrt(DDD)/(2*3^2));
% Exponential kernel
% DD=auto(GC.inverted.D)*auto(GC.inverted.D)';
% Gauss_kernel=(DD +1).^2;
% Sigmoid kernel
% DD=auto(GC.inverted.D)*auto(GC.inverted.D)';
% Gauss_kernel=tanh(10*DD -1);


function flag=validateCoeffcients(x)
c1=x(1,:);
c2=x(2,:);
c3=x(3,:);
coefficients1 = [1; -c1];
coefficients2 = [1; -c2];
coefficients3 = [1; -c3];
if max(abs(roots(coefficients1))) < 1 && max(abs(roots(coefficients2))) < 1 && max(abs(roots(coefficients3))) < 1
    flag = true;
else
    flag= false;
    warning('coefficients not stable')
end

end
