function [Xb,trace]=IDE_VTRGMM(D, DL, NP, gen_max, X, Y, Pi, beta, m, v, W, alpha0, m0, beta0, W0, v0,CR,F)

trace=zeros(gen_max,2);
x=floor(rand(NP,D)*DL);


cost=zeros(1,NP);
cost(1)=fitness_GMM(x(1,:), X, Y, DL, Pi, beta, m, v, W, alpha0, m0, beta0, v0, W0);
Pb=cost(1);
Xb=x(1,:);
for i=2:NP
    cost(i)=fitness_GMM(x(i,:), X, Y, DL, Pi, beta, m, v, W, alpha0, m0, beta0, v0, W0);
    if(cost(i)<=Pb)
        Pb=cost(i);
        Xb=x(i,:);
    end
end
trace(1,1)=1;
trace(1,2)=Pb;

count=1;
for j=1:gen_max
    
    for i=1:NP
        permuteFlag=true;
        while permuteFlag
            tempPerm=randperm(NP,3);
            a=tempPerm(1);
            b=tempPerm(2);
            c=tempPerm(3);
            if sum(tempPerm==i)==0
                permuteFlag=false;
            else
                permuteFlag=true;
            end
        end
        jrand=floor(rand*D+1);
        r=rand(D,1);
        trial_1=x(i,:);
        trial_2=x(i,:);
        for k=1:D
            crossoverFlag=r(k)<=CR;
            jFlag=jrand==k;
            if(crossoverFlag||jFlag)
                
                crossover=x(c,k)+F*(x(a,k)-x(b,k));
                if crossover <0
                    crossover=0;
                elseif crossover > DL
                    crossover=DL;
                end
                trial_1(k)=floor(crossover);
                trial_2(k)=ceil(crossover);
            end
        end
        score_1=fitness_GMM(trial_1, X, Y, DL, Pi, beta, m, v, W, alpha0, m0, beta0, v0, W0);
        score_2=fitness_GMM(trial_2, X, Y, DL, Pi, beta, m, v, W, alpha0, m0, beta0, v0, W0);
        
        if score_1<score_2
            trial = trial_1; 
            score=score_1;
        elseif score_2<score_1
            trial = trial_2;
            score=score_2;
        else
            trial = trial_1; 
            score=score_1;
        end
       
        if(score<=cost(i))
            x(i,1:D)=trial(1:D);
            cost(i)=score;
        end
        if cost(i)<=Pb
            Pb=cost(i);
            Xb(1:D)=x(i,1:D);
        end
    end
    
    count=count+1;
    trace(count,1)=count;
    trace(count,2)=Pb;
end

end
