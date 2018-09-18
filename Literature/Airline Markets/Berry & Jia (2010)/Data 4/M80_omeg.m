function omeg=M80_omeg(dM,p_mc,sH)
% dM has the following fields: dM.expmv, dM.XBxi
% dM.Midx, dM.MidxL, dM.Mid, dM.Cidx, dM.Sidx, dM.Sid, dM.theta, dM.s_jm

% p=c-q/(partial q/ partial p)
% dq/dp = (alpha/lambda) * s_j * [1 - s_jgrp +lambda*s_j/[sum_k_exp(delta/lambda)]^lambda]

%Parameters
alpha=dM.alpha;
lambda=dM.lambda;
gamma=dM.gamma;
ntype=dM.ntype;

%%% Markup for multi-product firms
omeg=zeros(dM.nobs,1);
for i=1:dM.nC
    id1=dM.CidxL(i)+1;
    id2=dM.CidxL(i+1);
    tn=id2-id1+1;
    
    tdqdp=zeros(tn,tn);
    for k=1:ntype
        temp1=sH.J(id1:id2,k);
        temp2=sH.jgrp(id1:id2,k)';
        tM=(lambda*sH.oth(id1,k)-1)* (temp1*temp2) + diag(temp1);
        tdqdp=tdqdp+ alpha(k)*gamma(k)*tM;
    end;
    tdqdp=tdqdp/lambda;
    % dqdp{i}=tdqdp;
    tq=dM.s_jm(id1:id2);
    omeg(id1:id2)=p_mc(id1:id2) + tdqdp \ tq;
end; 

