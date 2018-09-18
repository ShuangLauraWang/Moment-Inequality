diary P130_est.out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameters for year 1999
% Panle Jia 5/20/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short g

% Parameters to be adjusted depending on the model
ntype=2;            
ncoef=3;
yr=1999;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part One: Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load P130_MkDist data   %ProdID,dist,SmDist,MdDist,LgDist
data=sortrows(data,1);  %sort by ProdID
SmDist=data(:,3); MdDist=data(:,4); LgDist=data(:,5);
SmMdDist=SmDist+MdDist;

%%% 'data' in P100_db1b.mat contains the following fields:
% MkID CarrID ProdSh BinFare NoConn 
% dept Hub Hubdest Dist NumDest 
% NdestEnd FLLAS NoRoute HubConn Slot 
% NoCarr NoLcc AvgPop Phat PassOrg
% PassDest CarrOrgSh CarrDestSh CityOrg CityDest
% Avr_Direct Avr_Hub Avr_Hubdest deptalt P25
% P75	SeatT CommT Del_15 Del_30
% DumLcc, DirHat, ConnHat;
load P100_db1b           
nJ=nJ(:,2);

data=sortrows(data,39);    %Sort by MkID, CarrID, Fare
s_jm=data(:,3);
nconn=data(:,5);
dept=data(:,29)/90;   %alternative dept 
dirIV=data(:,37)/90;
connIV=data(:,38)/90;
deptIV=(nconn==0).*dirIV+(nconn==2).*connIV;

fare25=data(:,30)/100;
fare75=data(:,31)/100;

%%% Carrier Dummies
CarrID=data(:,2);
AA= (CarrID==1);
CO= (CarrID==3);
DL= (CarrID==4);
HP= (CarrID==9);
NW= (CarrID==13);
TW= (CarrID==14);
UA= (CarrID==16);
US= (CarrID==17);
WN= (CarrID==18);
OT=(CarrID~=1 &CarrID~=3 &CarrID~=4 &CarrID~=9 &CarrID~=13 &CarrID~=14 &CarrID~=16 &CarrID~=17 &CarrID~=18); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part Two: IV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fare=data(:,4)/100;        %in hund 
dist=data(:,9)/1000;      %in thou of miles
dist2=dist.^2;
hub=data(:,7); hubdest=data(:,8);
ndest=data(:,10)/100;      %in hund of connected cities
ndestend=data(:,11)/100;
hubconn=data(:,14);
lcc=data(:,36); 
tour=data(:,12);           %tour is dummy for FL / LAS
slot=data(:,15); pop=data(:,18)/1000000; 
nrout=data(:,13); nocarr=data(:,16); 
adir=data(:,26); 

slotMC= (slot>0); 
hubMC= (hub+hubdest+hubconn)>0;

% Demand IV
IV1=[ones(nobs,1),dist,nconn,ndest,ndestend, hub,hubdest,lcc,...   
    tour,slot,pop,nrout,nocarr,...     
    hubconn,...
    dist.*ndest, dist.*ndestend, dist.*lcc,... 
    dist.*tour, dist.*pop, ...
    nconn.*ndest, nconn.*tour, nconn.*pop, ...
    nconn.*nrout,  ...
    ndest.*lcc,  ndest.*pop, ndest.*nrout,...
    lcc.*tour, lcc.*pop,...
    tour.*nrout, tour.*pop, slot.*nrout,  pop.*nrout,  ...
    adir.*dist, adir.*nconn,  adir.*ndest, adir.*hubdest,...
    dist.*slot, nconn.*slot, ndest.*slot,ndestend.*slot, slot.*pop, ...
    fare25, fare75,OT,CO,DL,NW,HP,TW,UA,US,WN,...
    ndestend.*nconn, ndestend.*ndest, ndestend.*hub,ndestend.*lcc,...
    ndestend.*tour,ndestend.*pop,ndestend.*nrout,ndestend.*adir,...
    deptIV, deptIV.*ndestend, deptIV.*dist];

% Supply IV
IV2=[ndest,ndestend,hub,hubdest,lcc,...   
    tour,slot,pop,nrout,nocarr,...
     dist.*ndest, dist.*ndestend, dist.*lcc,... 
    dist.*tour, dist.*pop, ...
     nconn.*ndest, nconn.*lcc, nconn.*tour, nconn.*pop, ...
    nconn.*nrout,  ...
    ndest.*lcc,  ndest.*pop, ndest.*nrout,...
    lcc.*tour, lcc.*pop,...
    tour.*nrout, tour.*pop, slot.*nrout,  pop.*nrout,...
    dist.*slot, nconn.*slot, ndest.*slot,ndestend.*slot, slot.*pop,...
    fare25, fare75,OT,CO,DL,NW,HP,TW,UA,US,WN,...
    ndestend.*nconn, ndestend.*ndest, ndestend.*hub,ndestend.*lcc,...
    ndestend.*tour,ndestend.*pop,ndestend.*nrout,ndestend.*adir,...
    deptIV, deptIV.*ndestend, deptIV.*dist,slotMC,hubMC,...
    ones(nobs,1).*SmMdDist,dist.*SmMdDist,nconn.*SmMdDist,...
    ones(nobs,1).*LgDist, dist.*LgDist, nconn.*LgDist];

invA1 = inv(IV1'*IV1);     
invA2 = inv(IV2'*IV2);    

corrIV1=corr(IV1(:,2:end));
for i=1:size(IV1,2)-1
    corrIV1(i,i:end)=0;
    corrIV1(i,1:i-1)=corrIV1(i,1:i-1).*(abs(corrIV1(i,1:i-1))>0.85);
end;
display('max corr')
max(abs(corrIV1))'
display('mean IV1')
mean(IV1)'

corrIV2=corr(IV2(:,2:end));
for i=1:size(IV2,2)-1
    corrIV2(i,i:end)=0;
    corrIV2(i,1:i-1)=corrIV2(i,1:i-1).*(abs(corrIV2(i,1:i-1))>0.85);
end;
display('max corr')
max(abs(corrIV2))'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part Three: first stage estimation
% Theta=[alpha1,layover1, alpha2,layover2, 
% alpha3,layover3,
% const, ndest, dept, dist,dist2, tour, slot,delay,
% AA,CO,DL,HP,NW,TW,UA,US,WN
% cost_const, cost_dist, cost_layover, hub,slot,seat, comm, AA,CO,DL,HP,NW,TW,UA,US,WN,
% lambda, gamma1, gamm2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coef similar to theta in P10_est32d
theta0=[-0.737656356210933;-0.542245738886397;-6.05107946208118;-0.0704985812292955;-0.339917044024353;-8.67715414904550;0.356044105497042;0.0492842335026210;0.260012038262847;-0.0440319040860364;0.336389405456297;-0.206892236068573;-0.0817900261964014;-0.211578868479930;-0.0993915445518525;-0.129885770806569;-0.111480619881294;-0.122903775662614;0.163399122688848;-0.166531914853298;0.0627915208737847;0.856435904866960;0.254125389332529;-0.0661355521434858;1.38546150029252;0.0836097003418259;-0.0907841121922663;0.000235074924175214;0.0580937620440613;-0.0256153684790905;-0.0379930498525333;-0.0943462258668135;-0.179866058260164;-0.0203217348017060;0.0116382566913122;-0.0323263904034479;-0.0831693009133002;-0.126443961069628;0.785763805918461;0.624002455242966;];
lb=[-3;-3;-9;-3;-3;-11;
    -ones(6,1);
    -ones(9,1);
    -ones(8,1);
    -ones(9,1);
    0.1; 0.01; ];
ub=[-0.001;-.001;-2; -.001;-.001;-2;
    1; 0.5; 1; 0.1; 1; 1;    %restrict coef of dist2 and dept
    ones(9,1);
    2;1;1;4;1;1;1;1;
    ones(9,1);
    0.95;0.95; ];
% A and b are linear constraints
A=[]; b=[];

% use fmincon to refine the search
options=optimset('Display','iter','MaxIter',1000,'MaxFunEvals',1000,'GradObj','on',...
    'DiffMinChange',1e-6,'DerivativeCheck','on');
XMat=[fare,nconn, ones(nobs,1),ndest,dept, dist,dist2, tour,slot, OT,CO,DL,HP,NW,TW,UA,US,WN,...
    ones(nobs,1).*SmMdDist, dist.*SmMdDist, nconn.*SmMdDist,...
    ones(nobs,1).*LgDist, dist.*LgDist, nconn.*LgDist, hubMC,slotMC];
dM.MidxL=[0;Midx]; dM.Mid=Mid; 
dM.CidxL=[0;Cidx]; 
%%%%%******** DM.CIDX2 IMPORTANT *******%%%%%%%%%%
dM.Cidx2=nJC(:,1);      %Cidx2 is the MkID
dM.s_jm=s_jm; 
dM.xiold=zeros(nobs,1); 
dM.tol=1e-12;    %tighter criterion for the xi inversion 
dM.nM=nM; dM.nC=nC; dM.nobs=nobs; 
dM.tolL=exp(-dM.tol); dM.tolH=exp(dM.tol);
dM.ntype=ntype; dM.ncoef=ncoef;      
dM.modspec=68; dM.yr=yr;      % cost coef differ for Sm/Md and long routes
VM.iv1=IV1; VM.inv1=invA1; VM.iv2=IV2; VM.inv2=invA2;

clear CarrID AA CO DL HP NW TW UA US WN
clear fare* nconn time dist* hubdest ndestend hubconn tour slot nrout nocarr passorg passdest 
clear adir X Cidx Mid Midx 
clear dept_dir deptC s_jm temp* xiold 
clear IV* deptIV FareIV connIV dirIV dept dirdept dist slotMC

tic
[theta, fval,exitflag,output,laglamb,fgrad] = fmincon(@M130_gmmGredMM,theta0,...
    A,b,[],[],lb,ub,[],options,XMat,VM,dM,0)
display('it takes fmincon ... minutes')
comp_t = toc/60

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance of the parameter estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VarT xi omeg] = M130_varMM(theta,XMat,VM,dM,0);
disp('mean/std of xi')
[mean(xi) std(xi)]
disp('mean/std of omeg')
[mean(omeg) std(omeg)]

disp('theta0, theta, std')
[theta0 theta sqrt(diag(VarT))]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the optimal weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ng1=size(VM.iv1,2);
ng2=size(VM.iv2,2);
V1=zeros(ng1+ng2,ng1+ng2);
for i=1:nobs
    t1=VM.iv1(i,:)';
    t2=VM.iv2(i,:)';
    t=[xi(i)*t1;omeg(i)*t2];
    V1=V1+t*t';
end;
V1=inv(V1/nobs);
OptWt=V1;
save P130_est theta* fval* exitflag* output* VarT* OptWt


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part Six: second stage estimation
% Theta=[const1, alpha1, layover1, ndest1, swjb1, time, dept, lcc, tour,
% const2, alpha2,
% layover2, ndest2, swjb2, cost_const, cost_dist, cost_layover, cost_lcc, spk_cost,spk_dist, lambda, gamma]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VM.OptWt=OptWt;
tic
[theta2, fval2,exitflag2,output2,laglamb2,fgrad2] = fmincon(@M130_gmmGredMM,theta,...
    A,b,[],[],lb,ub,[],options,XMat,VM,dM,1)
display('it takes fmincon ... minutes')
comp_t2 = toc/60

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance of the parameter estimates
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dM.modspec=50;
[VarT2 xi2 omeg2] = M150_varMkDep(theta2,XMat,VM,dM);
disp('mean/std of xi2')
[mean(xi2) std(xi2)]
disp('mean/std of omeg2')
[mean(omeg2) std(omeg2)]
disp('theta, theta2, std2')
[theta theta2 sqrt(diag(VarT2))]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save P130_est theta* fval* exitflag* out* VarT* OptWt xi2 omeg2
diary off

