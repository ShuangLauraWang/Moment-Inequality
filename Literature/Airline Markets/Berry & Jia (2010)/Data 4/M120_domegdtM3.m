function domegdt=M120_domegdtM3(XMat,dM,sH,Kx,dxidt,dS,mkup)
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

%%%%% Parameters
alpha=dM.alpha; lambda=dM.lambda; gamma=dM.gamma;
ntype=dM.ntype; ncoef=dM.ncoef;
nobs=dM.nobs;

%%%%% Calculate dfdb
dM.Mid=[]; dM.MidxL=[];
dM.xiold=[]; dM.s_jm=[];
dM.expmv=[]; 
sH.spk=[];

sH.jgrp=sH.jgrp';
domegdt=zeros(nobs,(ncoef+1)*ntype);
for i=1:dM.nC
    id1=dM.CidxL(i)+1;
    id2=dM.CidxL(i+1);
    tn=id2-id1+1;   %number of products for the same carrier
    
    %dqdp=zeros(tn,tn);
    dqdb=zeros(tn,tn*((ncoef+1)*ntype));
    tmkup=mkup(id1:id2);
    id = dM.Cidx2(i)-1;
    
    tqtg=zeros(tn,ntype*tn);
    tdqdp=zeros(tn,tn);
    for j=1:ntype
        ts0=sH.oth(id1,j);
        lams0=lambda*ts0;

        ts1=sH.J(id1:id2,j);
        ts2=sH.jgrp(j,id1:id2);
        ts1s2=ts1*ts2;
        
        tmp=alpha(j)/lambda*( (lams0-1)*ts1s2+diag(ts1));
        tqtg(:,(j-1)*tn+1:j*tn)=tmp;
        tdqdp=tdqdp+gamma(j)*tmp;
        
        kx2=Kx(id1,4*j-3);  %fare
        kx3=Kx(id1,4*j-2);  %nconn
        kx3b=Kx(id1,4*j-1); %tour
        kx1=Kx(id1,4*j);    %const
        kx4=Kx(id1,4*ntype+j);
        
        tx2=XMat(id1:id2,1);    %fare
        tx3=XMat(id1:id2,2);    %nconn
        tx3b=XMat(id1:id2,3);   %tour
        tx1=XMat(id1:id2,4);    %const

        txM2=tx2(:,ones(1,tn));
        txM2=(txM2+txM2')/lambda;
        txM3=tx3(:,ones(1,tn));
        txM3=(txM3+txM3')/lambda;
        txM3b=tx3b(:,ones(1,tn));
        txM3b=(txM3b+txM3b')/lambda;
        txM1=tx1(:,ones(1,tn));
        txM1=(txM1+txM1')/lambda;

        for k = 1: ((ncoef+1)*ntype)
            tdxidt=dxidt(id1:id2,k);
            tdel=dS.J(ntype*id+j,k);
            tdel2=-dS.jgrp(ntype*id+j,k)/lambda;
            tdel3=dS.oth(ntype*id+j,k);
            tdel_diag=tdel + 1/lambda*tdxidt;

            txi=tdxidt(:,ones(1,tn));
            txi=(txi+txi')/lambda;
            tdelM=tdel+tdel2+txi;

            %partial q1/partial const1,beta1
            if k==4*j-3     %partial q/partial alpha
                c=(lams0-1)*( txM2 + kx2*(lams0-2) + tdelM );
                cc=lambda*(1-ts0)*kx2+tdel3;
                tmpdiag=ts1.*(tx2/lambda + kx2*(lams0-1) + tdel_diag);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k) + gamma(j)*alpha(j)/lambda* (ts1s2.*(c - lams0*cc) + diag(tmpdiag) )+...
                    gamma(j)/lambda* ( (lams0-1)*ts1s2 + diag(ts1) );

            elseif k==4*j-2    %partial q/partial nconn
                c=(lams0-1)*( txM3 + kx3*(lams0-2) + tdelM );
                cc=lambda*(1-ts0)*kx3+tdel3;
                tmpdiag=ts1.*(tx3/lambda + kx3*(lams0-1) + tdel_diag);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k)+gamma(j)*alpha(j)/lambda* (ts1s2.*(c - lams0*cc) + diag(tmpdiag) );
            elseif k==4*j-1    %partial q/partial tour
                c=(lams0-1)*( txM3b + kx3b*(lams0-2) + tdelM );
                cc=lambda*(1-ts0)*kx3b+tdel3;
                tmpdiag=ts1.*(tx3b/lambda + kx3b*(lams0-1) + tdel_diag);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k)+gamma(j)*alpha(j)/lambda* (ts1s2.*(c - lams0*cc) + diag(tmpdiag) );

            elseif k==4*j    %partial q/partial const
                c=(lams0-1)*( txM1 + kx1*(lams0-2) + tdelM );
                cc=lambda*(1-ts0)*kx1+tdel3;
                tmpdiag=ts1.*(tx1/lambda + kx1*(lams0-1) + tdel_diag);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k)+gamma(j)*alpha(j)/lambda* (ts1s2.*(c - lams0*cc) + diag(tmpdiag) );

            elseif k==(ncoef*ntype+1) %partial q1/partial lambda
                tmv=-dM.XBxi(id1:id2,j)/(lambda^2);
                c0=ts0*sH.lnsum(id1,j)+kx4*(lams0-1);
                tdiag=ts1.*( tmv+c0+tdel_diag-1/lambda ) /lambda;

                tmv=tmv(:,ones(1,tn));
                tmv=tmv+tmv';
                c1=(ts0-1/lambda)*(tmv + c0-kx4 + tdelM);
                c2=1/(lambda^2)-ts0*( (1-ts0)*(sH.lnsum(id1,j)+lambda*kx4) + tdel3);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k)+gamma(j)*alpha(j)*( ts1s2.* (c1 + c2) + diag(tdiag) );
            else 
                c=(lams0-1)* tdelM ;
                tmpdiag=ts1.*(tdel_diag);
                dqdb(:,tn*(k-1)+1:tn*k)=dqdb(:,tn*(k-1)+1:tn*k)+gamma(j)*alpha(j)/lambda* (ts1s2.* (c - lams0*tdel3) + diag(tmpdiag) );
            end;
        end;
        
    end;
    % partial q / partial gamma
    tmp=tqtg(:,(ntype-1)*tn+1:end); tmp=repmat(tmp,1,ntype-1);
    tqtgM=tqtg(:,1:(ntype-1)*tn)-tmp;
    dqdb(:,(ncoef*ntype+1)*tn+1:end)=dqdb(:,(ncoef*ntype+1)*tn+1:end)+tqtgM;
    
    %%% dfdb is 0 for betas common to both types
    dfdb=zeros(tn,(ncoef+1)*ntype);
    dfdb(:)=tmkup' * dqdb;
    domegdt(id1:id2,:)= tdqdp \ dfdb;
end;

