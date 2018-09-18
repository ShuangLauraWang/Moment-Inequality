function [VarT xi omeg] = M130_varMM(theta,XMat,VM,dM,stg2)
% The file should accommodate changes in specifications for both years

% XMat=[fare,nconn,ones(nobs,1), ndest,swjb,dept, dist, dist2,lcc,tour,slot, dist,nconn,(hub/slot)];
% VM.iv1=IV1; VM.inv1=invA1; VM.iv2=IV2; VM.inv2=invA2;
% VM.iv3=IV3; VM.inv3=invA3; 

% dM has the following fields: 
% dM.MidxL, dM.Mid, dM.CidxL,
% dM.SidxL, dM.Sid,
% dM.s_jm, dM.xiold, dM.tol, 
% dM.nM, dM.nC, dM.nS, dM.nobs,
% dM.tol, dM.tolH, dM.tolL, which should be exp(+/-1e-7)

% Use the fact that: dxi/dt=-X for theta6-theta9 (the beta common to both types),
% dxi/dt=0 for theta15-18 (coef for product mc), 
% domeg/dt=0 for theta6-9, -X for theta15-18, 
% dS1/dS2/dxidt,domegdt has 14 columns: theta1-6, theta11-16, lambda, gamma
% stg2=0: first stage estimation; stg2=1: second stage estimation

modspec=dM.modspec; yr=dM.yr;
ntype=dM.ntype; ncoef=dM.ncoef; nobs=dM.nobs;
dM.lambda=theta(end-ntype+1); 
tmp=theta(end-ntype+2:end); dM.gamma=[tmp', 1-sum(tmp)];
[dM.XB dM.alpha p_mc]=M130_xbMM(XMat, theta,modspec,yr);

%%% Error from the demand equation: xi
[xi dM.expmv] = M80_xi(dM);
dM.XBxi=dM.XB+xi(:,ones(1,ntype));
dM.XB=[]; dM.xiold=[];

%%% Error from optimal prices: omeg
sH=M75_sh(dM);
omeg=M80_omeg(dM,p_mc,sH);
mkup=p_mc-omeg;

Kx=M80_KxM(dM,XMat(:,1:ncoef));
[dxidt,dS]=M80_dxidtM(XMat(:,1:ncoef),dM,sH,Kx);
if ncoef==3
    domegdt=M80_domegdtM3(XMat(:,1:ncoef),dM,sH,Kx,dxidt,dS,mkup);
elseif ncoef==4
    domegdt=M120_domegdtM3(XMat(:,1:ncoef),dM,sH,Kx,dxidt,dS,mkup);
end;            
clear sH Kx dS dM mkup

if yr==2006
    switch modspec
        case 22     % 2 types, 3 coef (include a type-specific const), 1 dept, 
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,13), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[3,18:21,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 68     % cost coef differs for short-medium/long haul routes
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,16), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 69     % 2 sets of cost coef, 2 types, 4 coef (include a type-specific const, tour=0 for busi), 1 dept, 
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):17), zeros(nobs,16), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 70     % 2 sets of cost coef, 2 types, 4 coef (include a type-specific const, tour), 1 dept, 
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,16), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 72     % 2 sets of cost coef, 2 types, 4 coef (include a type-specific const, tour=0 for busi), 1 dept, 
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):16), zeros(nobs,15), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,16-ncoef),-XMat(:,[17:24,10:16]),domegdt(:,ncoef*ntype+1:end)];
        case 73     % delay in demand, seat/delay in supply
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):18), zeros(nobs,18), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,18-ncoef),-XMat(:,[19:28,11:18]),domegdt(:,ncoef*ntype+1:end)];
        case 74     % 25 airports
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):17), zeros(nobs,16),-XMat(:,26:50), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),zeros(nobs,25),domegdt(:,ncoef*ntype+1:end)];
        case 75     % Mk>3k,tour=0 for busi
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):17), zeros(nobs,13), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,17-ncoef),-XMat(:,[18:22,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 76     % MK>3k, no B6/SW Entry, tour=0 for busi
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):16), zeros(nobs,12), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,16-ncoef),-XMat(:,[17:21,10:16]),domegdt(:,ncoef*ntype+1:end)];
        case 77     % Mk>3k, 2 tour coef
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,13), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:22,10:17]),domegdt(:,ncoef*ntype+1:end)];
        case 78     % Mk>3k, 1 tour coef, 1 dist
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):16), zeros(nobs,13), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,16-ncoef),-XMat(:,[17:21,9:16]),domegdt(:,ncoef*ntype+1:end)];
    end;
elseif yr==1999
    switch modspec
        case 22     % 2 types, 3 coef (including a type specific const), 1 dept, 
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,14), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[3,19:22,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 68     % cost coef differ for short_medium/long haul routes
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 69     % 2 sets of cost coef, 2 types, 4 coef (include a type-specific const, tour=0 for busi), 1 dept, 
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):18), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 70     % 2 sets of cost coef, 2 types, 4 coef (include a type-specific const, tour), 1 dept, 
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 73     % delay in demand, delay/seat in supply
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):19), zeros(nobs,19), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,19-ncoef),-XMat(:,[20:29,11:19]),domegdt(:,ncoef*ntype+1:end)];
        case 74     % 25 airports
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):18), zeros(nobs,17), -XMat(:,27:51), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),zeros(nobs,25),domegdt(:,ncoef*ntype+1:end)];
        case 75     % Mk longer than 3k,tour=0 for busi
            dxidt=[dxidt(:,[1:6,8]),-XMat(:,(ncoef+1):18), zeros(nobs,14), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,[1:6,8]),zeros(nobs,18-ncoef),-XMat(:,[19:23,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 77     % Mk>3k, 2 tour coef
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,14), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:23,10:18]),domegdt(:,ncoef*ntype+1:end)];
        case 78     % Mk>3k, 1 tour coef, 1 dist
            dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,14), dxidt(:,ncoef*ntype+1:end)];
            domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:22,9:17]),domegdt(:,ncoef*ntype+1:end)];
    end;
end;            

% Variance matrix of theta
ng1=size(VM.iv1,2); ng2=size(VM.iv2,2);
Gamma=[VM.iv1'*dxidt/nobs; VM.iv2'*domegdt/nobs];
if stg2==0
    V1=zeros(ng1+ng2,ng1+ng2);
    for i=1:nobs
        t1=VM.iv1(i,:);
        t2=VM.iv2(i,:);
        t=[xi(i)*t1, omeg(i)*t2];
        V1=V1+t'*t;
    end;
    V1=V1/nobs;

    %Cphi is used in the first stage to mimic optimal Weight matrix
    V=V1/nobs;
    Wgt=blkdiag(VM.inv1, VM.inv2);
    Beta=Gamma'*Wgt*Gamma;
    Sigma=Gamma'*Wgt*V*Wgt*Gamma;
    VarT=inv(Beta)*Sigma*inv(Beta);
elseif stg2==1
    Beta=Gamma'*VM.OptWt*Gamma;
    VarT=inv(Beta)/nobs;
end;

