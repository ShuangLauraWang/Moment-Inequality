function [BetaPoint, fval, ErrorFlag]=f_fminconsingle(ZJ,WJ,lb,ub,options,Beta_Start);

    ErrorFlag = 1;
    NumBetaStart = size(Beta_Start,2);
    fval = 10^20;

    for n=1:NumBetaStart
        [BetaTemp, fvalTemp] = fminsearch(@f_singlepointmoment,Beta_Start(:,n),options,ZJ,WJ);
        if (fvalTemp < fval)
            fval = fvalTemp;
            BetaPoint = BetaTemp;
            ErrorFlag = 0;
        end;
    end;
    if (ErrorFlag)
        BetaPoint = zeros(size(Beta_Start(:,1)),1);
    end;