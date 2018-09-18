function [VarT xi omeg] = M150_varMkDep(theta,XMat,VM,dM)
% Allowing arbitrary dependence among products in the same market
% dM.CidxL necessary to determine the group of products with dependence

modspec=dM.modspec; yr=dM.yr;
ntype=dM.ntype; ncoef=dM.ncoef; nobs=dM.nobs;
dM.lambda=theta(end-ntype+1); 
tmp=theta(end-ntype+2:end); dM.gamma=[tmp', 1-sum(tmp)];
[dM.XB dM.alpha p_mc]=M150_xbMM(XMat, theta,modspec,yr);

%%% Error from the demand equation: xi
[xi dM.expmv] = M80_xi(dM);
dM.XBxi=dM.XB+xi(:,ones(1,ntype));
dM.XB=[]; dM.xiold=[];

%%% Error from optimal prices: omeg
sH=M75_sh(dM);
omeg=M80_omeg(dM,p_mc,sH);
mkup=p_mc-omeg;

% Obtain the optimal weight
ng1=size(VM.iv1,2);
ng2=size(VM.iv2,2);
V1=zeros(ng1+ng2,ng1+ng2);
for i=1:dM.nM
    id1=dM.MidxL(i)+1; id2=dM.MidxL(i+1);
    t1=VM.iv1(id1:id2,:);
    t2=VM.iv2(id1:id2,:);
    
    txi=xi(id1:id2); txi=repmat(txi,1,ng1);
    tomeg=omeg(id1:id2); tomeg=repmat(tomeg,1,ng2);
    tg=[txi.*t1, tomeg.*t2];
    tg=sum(tg,1);       %important to use sum(tg,1) instead of sum(tg)
    V1=V1+tg'*tg;
end;
OptWt=inv(V1/nobs);
VM.OptWt=OptWt;

Kx=M80_KxM(dM,XMat(:,1:ncoef));
[dxidt,dS]=M80_dxidtM(XMat(:,1:ncoef),dM,sH,Kx);
domegdt=M80_domegdtM3(XMat(:,1:ncoef),dM,sH,Kx,dxidt,dS,mkup);
clear sH Kx dS dM mkup

    if yr==2006
        switch modspec
            case 30     % no LCC
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):16), zeros(nobs,15), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,16-ncoef),-XMat(:,[17:24,10:16]),domegdt(:,ncoef*ntype+1:end)];
            case 44     % delay in demand only
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,16), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,11:18]),domegdt(:,ncoef*ntype+1:end)];
            case 50     % combine airports; same as basecase
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,16), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),domegdt(:,ncoef*ntype+1:end)];
            case 80     % 25 airport dummies
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,16),-XMat(:,26:50), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:25,10:17]),zeros(nobs,25),domegdt(:,ncoef*ntype+1:end)];
            case 90     %Mk>1.5k, 1 set of cost, 1 dist, tour same for both types 
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):16), zeros(nobs,13), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,16-ncoef),-XMat(:,[17:21,9:16]),domegdt(:,ncoef*ntype+1:end)];
        end;
    elseif yr==1999
        switch modspec
            case 30     % no LCC
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),domegdt(:,ncoef*ntype+1:end)];
            case 44     % delay in demand only
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):19), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,19-ncoef),-XMat(:,[20:27,11:19]),domegdt(:,ncoef*ntype+1:end)];
            case 50     % combine airports; same as basecase
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,17), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),domegdt(:,ncoef*ntype+1:end)];
            case 80     % 25 airport dummies
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):18), zeros(nobs,17),-XMat(:,27:51), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,18-ncoef),-XMat(:,[19:26,10:18]),zeros(nobs,25),domegdt(:,ncoef*ntype+1:end)];
            case 90     %Mk>1.5k, 1 set of cost, 1 dist, tour same for both types 
                dxidt=[dxidt(:,1:ncoef*ntype),-XMat(:,(ncoef+1):17), zeros(nobs,14), dxidt(:,ncoef*ntype+1:end)];
                domegdt=[domegdt(:,1:ncoef*ntype),zeros(nobs,17-ncoef),-XMat(:,[18:22,9:17]),domegdt(:,ncoef*ntype+1:end)];
        end;
    end;            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Variance matrix of theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The gradient and covariance matrix
Gamma=[VM.iv1'*dxidt/nobs; VM.iv2'*domegdt/nobs];
Beta=Gamma'*VM.OptWt*Gamma;
VarT=inv(Beta)/nobs;

