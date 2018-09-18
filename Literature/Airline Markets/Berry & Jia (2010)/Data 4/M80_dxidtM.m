function [dxidt,dS]=M80_dxidtM(XMat,dM,sH,Kx)
% dxidt=-X for the betas that are the same across the two groups

% XMat=[ones(nobs,1),fare, noconn,hub, swjb];
% dM has the following fields: % dM.Midx, dM.MidxL, dM.Mid, dM.CidxL,
% dM.Sidx, dM.Sid, dM.s_jm, dM.xiold, dM.tol, dM.nM, dM.nC, dM.nobs,
% dM.X1/X2, dM.expmv1/dM.expmv2
% dM.tol, dM.tolH, dM.tolL

% sH1 has the following fields: sH1.J, sH1.jgrp, sH1.oth, sH1.lnsum
% sH2 has the following fields: sH2.J, sH2.jgrp, sH2.oth, sH2.lnsum
% Kx1: const, alpha1, b_layover, b_hub, b_swjb, b_dept, lamb
% Kx2: const, alpha2, b_layover, b_hub, b_swjb, b_dept, lamb

% dS1 has fields: dS1.jgrp=ds1grpdxidt, dS1.J=ds1dxidt, dS1.oth=ds10dxidt;
% dS2 has fields: dS2.jgrp, dS2.J, dS2.oth;

% Use the fact that: dxi/dt=-X for theta7-theta10 (the beta common to both types),
% dxi/dt=0 for theta17-20 (coef for product mc), 
% domeg/dt=0 for theta7-10, -X for theta17-20, 
% dS1/dS2/dxidt,domegdt has 14 columns: theta1-6, theta11-16, lambda, gamma

% Remove irrelevant fields to save memory
dM.Midx=[]; dM.CidxL=[]; dM.Cidx2=[]; dM.Sidx=[]; dM.Sid=[]; dM.s_jm=[]; dM.xiold=[];
dM.expmv=[]; 
lambda=dM.lambda; gamma=dM.gamma;
ntype=dM.ntype; ncoef=dM.ncoef;
nobs=dM.nobs; modspec=dM.modspec;

%%% partial_xi/partial_theta
idx1=kron(1:ntype,ones(1,ncoef));
idx2=repmat(1:ncoef,1,ntype);
gammaM=gamma(ones(1,nobs),idx1);

lams0M=lambda*sH.oth(:,idx1);
dsdb=gammaM .* sH.J(:,idx1) .* ( XMat(:,idx2)/lambda + Kx(:,1:ncoef*ntype).*(lams0M-1) );
clear gammaM lams0M
dsdlam=sH.J .* (-dM.XBxi/(lambda^2) + sH.oth.*sH.lnsum + Kx(:,(ncoef*ntype)+1:(ncoef+1)*ntype).*(lambda*sH.oth-1));
dsdlam=M80_GamSum(dsdlam,gamma,ntype);
dsdgam=sH.J(:,1:ntype-1) - sH.J(:,ones(1,ntype-1)*ntype);

%%% dsdt: (ncoef+1)*ntype columns
dsdt=[dsdb, dsdlam, dsdgam];
clear dsdb dsdlam dsdgam XMat Kx 

dxidt=zeros(dM.nobs,(ncoef+1)*ntype);
if modspec~=41
    dS.jgrp=zeros(dM.nM*ntype,(ncoef+1)*ntype);
    dS.J=zeros(dM.nM*ntype,(ncoef+1)*ntype);
    dS.oth=zeros(dM.nM*ntype,(ncoef+1)*ntype);
elseif modspec==41
    dS=[];
end;

sH.oth=sH.oth'; sH.jgrp=sH.jgrp';
for i=1:dM.nM
    id1=dM.MidxL(i)+1;
    id2=dM.MidxL(i+1);
    tn=id2-id1+1;
    
    s0gam=sH.oth(:,id1)-1/lambda; 
    s_jgrp=sH.jgrp(:,id1:id2);
    s_j=sH.J(id1:id2,:);

    temp=zeros(tn,tn);
    for k=1:ntype
        temp=temp+ gamma(k) * (s0gam(k)*(s_j(:,k)*s_jgrp(k,:)) + diag(s_j(:,k))/lambda);
    end;
    tdsdt=dsdt(id1:id2,:);
    tdxidt= - temp \ tdsdt;
    
    if modspec~=41
        dS.jgrp(ntype*(i-1)+1:ntype*i,:)=s_jgrp*tdxidt;
        dS.J(ntype*(i-1)+1:ntype*i,:)= s0gam(:,ones(1,(ncoef+1)*ntype)) .* dS.jgrp(ntype*(i-1)+1:ntype*i,:);
        dS.oth(ntype*(i-1)+1:ntype*i,:)= (s_j') *tdxidt;
    end;        
    dxidt(id1:id2,:)=tdxidt;
end;

clear sH dM
