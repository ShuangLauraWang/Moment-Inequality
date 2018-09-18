function sH=M75_sh(dM)

% dM has the following fields: 
% dM.Midx, dM.MidxL, dM.Mid, dM.CidxL, dM.Sidx, dM.Sid,
% dM.s_jm, dM.xiold, dM.tol, dM.nM, dM.nC, dM.nS, dM.nobs,
% dM.tol, dM.tolH, dM.tolL

lamb=dM.lambda;
%ntype=dM.ntype;

tsum=M40_MkSum(dM.nM,dM.ntype,dM.nobs,dM.MidxL,dM.expmv);
if min(min(tsum))==0
    error('sum of shares = 0')
end;
grpsh=(tsum.^lamb)./(1+tsum.^lamb);
grpsh=grpsh(dM.Mid,:);

tsum = tsum(dM.Mid,:);
sH.lnsum=log(tsum);
sH.jgrp = dM.expmv./tsum;
sH.oth=1-grpsh;
sH.J = sH.jgrp.*grpsh;

clear dM grpsh tsum

