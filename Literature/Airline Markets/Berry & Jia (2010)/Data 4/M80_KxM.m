function Kx=M80_KxM(dM,XMat)

% dM has the following fields: dM.expmv1, dM.expmv2, dM.X1=X*beta1+xi,
% dM.X2=X*beta2+xi
% dM.cdindex, dM.cdindexL, dM.cdid, dM.cdindexCarr, dM.theta, dM.s_jm
% XMat: ones, fare, layover,hub,swjb,dept
% Kx1: const, alpha1, b_layover, b_hub, b_swjb, b_dept, lamb
% Kx2: const, alpha2, b_layover, b_hub, b_swjb, b_dept, lamb

% Use the fact that: dxi/dt=-X for theta7-theta10 (the beta common to both types),
% dxi/dt=0 for theta17-20 (coef for product mc), 
% domeg/dt=0 for theta7-10, -X for theta17-20, 
% dS1/dS2/dxidt,domegdt has 14 columns: theta1-6, theta11-16, lambda, gamma

lamb=dM.lambda;
ntype=dM.ntype; ncoef=dM.ncoef;
nobs=dM.nobs;
const=ones(nobs,ntype);

temp=[repmat(XMat(:,1:ncoef),1,ntype),-dM.XBxi/lamb,const];
clear XMat 
dM.XBxi=[]; 

idx1=kron(1:ntype,ones(1,ncoef));
temp=dM.expmv(:,[idx1,1:ntype,1:ntype]).*temp;
dM.expmv=[];

Kx=M40_MkSum(dM.nM, (ncoef+2)*ntype, dM.nobs, dM.MidxL, temp);

temp=Kx(:,(ncoef+1)*ntype+1:end);
Kx=Kx(:,1:(ncoef+1)*ntype)./temp(:,[idx1,1:ntype])/lamb;
Kx=Kx(dM.Mid,:);

clear dM temp

