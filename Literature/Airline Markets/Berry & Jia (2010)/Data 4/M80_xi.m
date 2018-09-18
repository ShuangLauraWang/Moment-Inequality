function [xi, Expmv] = M80_xi(dM)
% This function computes the mean utility level
% XMat: ones, fare, layover,hub,time,dept,lcc
% theta: a column vector

lambda=dM.lambda;
gamma=dM.gamma;
ntype=dM.ntype;

%%%%%%% INVERT FOR THE MEAN UTILITY GIVEN THETA
xi=zeros(dM.nobs,1);
Expmv = zeros(dM.nobs,ntype);     

%%% Needs to guard against 0 or inf
expXB=exp(dM.XB/lambda);
expXB=min(expXB,realmax/(1e+10));   %cap expXB at 1.8E+298
expXB=max(expXB,realmin*(1e+10));   %floor expXB at 2.2E-298

exiold=exp(dM.xiold/lambda);
for m=1:dM.nM
    id1=dM.MidxL(m)+1;
    id2=dM.MidxL(m+1);
    tn=id2-id1+1;
    
    teX=expXB(id1:id2,:);
    ts_jm=dM.s_jm(id1:id2); texiold=exiold(id1:id2);
    
    normH=10; normL=0.1;    
    while normH > dM.tolH || normL <dM.tolL
        expmv=teX.*texiold(:,ones(1,ntype));
        sumV=sum(expmv);        %a row vector. Should check whether only 1 prod in a market
        sumV=(sumV.^(lambda-1))./(1+sumV.^lambda);
        sumV=sumV(ones(tn,1),:);
        
        s_j = expmv.*sumV;
        mkshr=M80_GamSum(s_j,gamma,ntype);

        texi = texiold .* (ts_jm./mkshr); 
        t = texi ./ texiold;
        
        if min(t)==0 || sum(isnan(t))>0 || sum(isinf(t))>0 || isreal(t)==0 
            save prob.mat
            error('xi is 0, NaN, Inf, or complex')
        end;
        
        normH = max(t);
        normL = min(t);
        texiold = texi;
    end
    xi(id1:id2)=log(texi)*lambda;
    Expmv(id1:id2,:)=teX.*texi(:,ones(1,ntype));
end
