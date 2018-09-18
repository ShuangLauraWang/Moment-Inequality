diary M130_est.out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate parameters for year 2006
% Panle Jia 5/21/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
format short g

% Parameters to be adjusted depending on the model
ntype=2;            
ncoef=3;
yr=2006;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part One: Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load M130_MkDist    %ProdID,MkDist,SmDist,MdDist,LgDist
data=sortrows(data,1);  %sort by ProdID
SmDist=data(:,3); MdDist=data(:,4); LgDist=data(:,5);
SmMdDist=SmDist+MdDist; clear SmDist MdDist

%%% 'data' in M100_db1b.mat contains the following fields:
% MkID CarrID ProdSh BinFare NoConn 
% dept Hub Hubdest Dist NumDest 
% NdestEnd FLLAS NoRoute HubConn Slot 
% NoCarr NoLcc AvgPop Phat PassOrg
% PassDest CarrOrgSh CarrDestSh CityOrg CityDest
% Avr_Direct Avr_Hub Avr_Hubdest deptalt P25
% P75	SeatT CommT Del_15 Del_30
% DumLcc DirHat ConnHat ProdID;
load M100_db1b           
nJ=nJ(:,2);

data=sortrows(data,39);    %Sort by MkID, CarrID, ProdID
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
B6= (CarrID==4);
CO= (CarrID==5);
DL= (CarrID==6);
NW= (CarrID==11);
UA= (CarrID==15);
US= (CarrID==16);
WN= (CarrID==17);
OT= (CarrID~=1 &CarrID~=4 &CarrID~=5 &CarrID~=6 &CarrID~=11 &CarrID~=15 &CarrID~=16 &CarrID~=17);

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
    fare25, fare75,OT,CO,DL,NW,UA,US,WN,...
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
    fare25, fare75,OT,CO,DL,NW,UA,US,WN,...
    ndestend.*nconn, ndestend.*ndest, ndestend.*hub,ndestend.*lcc,...
    ndestend.*tour,ndestend.*pop,ndestend.*nrout,ndestend.*adir,...
    deptIV, deptIV.*ndestend, deptIV.*dist,slotMC,hubMC,...
    ones(nobs,1).*SmMdDist, dist.*SmMdDist, nconn.*SmMdDist,...
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
% const, ndest, dept, dist,dist2, tour, slot,
% OT,B6,CO,DL,NW,UA,US,WN
% cost_const, cost_dist, cost_layover, cost_lcc, lambda, gamma1, gamm2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0=[-0.942682104142090;-0.604539496550158;-6.00805822204041;-0.0938823173238298;-0.554825006228195;-8.91769138059189;0.267240103856126;0.110578136185204;0.512019372128040;-0.0775887159476948;0.395400127008418;-0.199709503628035;0.0784365751453704;0.447290565282301;0.0412530767886149;-0.219831799288678;0.0955608695060965;0.0709122737035883;0.0752477586604947;0.108152680634227;1.00676949929228;0.176647335633860;0.0810351688753235;1.40774437986201;0.0518240930242996;0.0535865547079051;-0.0509479293358409;0.00802671336972054;-0.166159188421474;-0.223566715371866;-0.151439301997620;-0.140108254547626;-0.0404793462828491;-0.0390045914247194;-0.0646696491322503;-0.141064041976670;0.742563699963633;0.544665953084210;];
lb=[-3;-3;-10;-3;-3;-10;
    -2*ones(6,1);
    -ones(8,1);
    -ones(8,1);
    -ones(8,1);
    0.1; 0.01];
ub=[-0.001;-0.001; -2; -0.001;-0.001; -2; 
    1; 0.5; 1; 0.5; 1; 0.1;
    ones(8,1);
    3;1;1;4;1;1; 1;1;
    ones(8,1);
    0.95;0.95];
% A and b are linear constraints
A=[]; b=[];

% use fmincon to refine the search
options=optimset('Display','iter','MaxIter',1000,'MaxFunEvals',1000,'GradObj','on',...
    'DiffMinChange',1e-6,'DerivativeCheck','on');
XMat=[fare,nconn, ones(nobs,1),ndest, dept, dist,dist2, tour,slot, OT,B6,CO,DL,NW,UA,US,WN,...
    ones(nobs,1).*SmMdDist, dist.*SmMdDist, nconn.*SmMdDist,...
    ones(nobs,1).*LgDist, dist.*LgDist, nconn.*LgDist, hubMC,slotMC];
dM.MidxL=[0;Midx]; dM.Mid=Mid; 
dM.CidxL=[0;Cidx]; 
%%%%****** IMPORTANT: DM.CIDX2 ***********%%%%
dM.Cidx2=nJC(:,1);      %MkID for each Mk/Carr combination
dM.s_jm=s_jm; 
dM.xiold=zeros(nobs,1); 
dM.tol=1e-12;    %tighter criterion for xi inversion than fmincon (outside loop)
dM.nM=nM; dM.nC=nC; dM.nobs=nobs; 
dM.tolL=exp(-dM.tol); dM.tolH=exp(dM.tol);
dM.ntype=ntype; dM.ncoef=ncoef;    
dM.modspec=68; dM.yr=yr;      %diff cost coef for short-medium/long routes
VM.iv1=IV1; VM.inv1=invA1; VM.iv2=IV2; VM.inv2=invA2;

clear CarrID AA CO DL NW UA US WN OT B6
clear fare* nconn time dist* hubdest ndestend hubconn tour slot nrout nocarr passorg passdest 
clear adir data Cidx Mid Midx 
clear dept_dir deptC s_jm temp* xiold 
clear IV* deptIV FareIV connIV dirIV dept dist slotMC hubMC regjet SmMdDist LgDist
 
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
save M130_est theta* fval* exitflag* output* VarT* OptWt


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

save M130_est theta* fval* exitflag* output* lag* fgrad* VarT* OptWt out* xi2 omeg2
diary off




