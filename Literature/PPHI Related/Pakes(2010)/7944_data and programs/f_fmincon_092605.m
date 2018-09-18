function [BetaMin, BetaMax, ErrorFlag]=f_fmincon_092605(ZJ,WJ,lb,ub,options,Beta_Start);
    ErrorFlag = 1;
    NumBetaStart = size(Beta_Start,2);
    DimBeta = size(Beta_Start,1);
    BetaMin=ub;
    BetaMax=lb;
    BetaMinTemp=zeros(DimBeta,NumBetaStart,1);
    BetaMaxTemp=zeros(DimBeta,NumBetaStart,1);    
    for m=1:DimBeta
        for n=1:NumBetaStart
%            display(Beta_Start(:,n));
            [BetaMinTemp(:,1),fval,exitflag] = fmincon(@f_calcbeta,Beta_Start(:,n),-ZJ,-WJ,[],[],lb,ub,[],options,ZJ,WJ,m,1);
                if ((exitflag >0)&&(BetaMinTemp(m,1)<BetaMin(m,1))) BetaMin(m,1)=BetaMinTemp(m,1); ErrorFlag=0; end;
            [BetaMaxTemp(:,1),fval,exitflag] = fmincon(@f_calcbeta,Beta_Start(:,n),-ZJ,-WJ,[],[],lb,ub,[],options,ZJ,WJ,m,0);
                if ((exitflag >0)&&(BetaMaxTemp(m,1)>BetaMax(m,1))) BetaMax(m,1)=BetaMaxTemp(m,1); ErrorFlag=0; end;
        end;
    end;
