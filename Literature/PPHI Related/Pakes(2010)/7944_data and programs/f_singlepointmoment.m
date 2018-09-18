function moment=f_singlepointmoment(Beta, ZJ, WJ);
    moment = 0;
    X = ZJ*Beta-WJ;
    n = size(X,1);
    for i=1:n
        if X(i,1) < 0
            moment = moment - X(i,1);
        end;
    end;