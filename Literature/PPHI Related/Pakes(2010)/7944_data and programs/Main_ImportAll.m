clear;

importdata=1;
if importdata
    x10 = importdata ('DATA/data.txt');
    x20 = importdata ('DATA/data_AllNetworkStructures.txt');

    x1a = importdata ('DATA/data_100306.txt');
    x2a = importdata ('DATA/data_AllNetworkStructures_100306.txt');
    x1b = importdata ('DATA/data_100806.txt');
    x2b = importdata ('DATA/data_AllNetworkStructures_100806.txt');
    x1c = importdata ('DATA/data_101206.txt');
    x2c = importdata ('DATA/data_AllNetworkStructures_101206.txt');
    x1d = importdata ('DATA/data_102206.txt');
    x2d = importdata ('DATA/data_AllNetworkStructures_102206.txt');
    x1e = importdata ('DATA/data_102206_server.txt');
    x2e = importdata ('DATA/data_AllNetworkStructures_102206_server.txt');
    x1g = importdata ('DATA/data_111306_server.txt');
    x2g = importdata ('DATA/data_AllNetworkStructures_111306_server.txt');
    x1h = importdata ('DATA/data_111406.txt');
    x2h = importdata ('DATA/data_AllNetworkStructures_111406.txt');
    x1i = importdata ('DATA/data_111506.txt');
    x2i = importdata ('DATA/data_AllNetworkStructures_111506.txt');
    x1j = importdata ('DATA/data_010707.txt');
    x2j = importdata ('DATA/data_AllNetworkStructures_010707.txt');
    x1k = importdata ('DATA/data_013007.txt');
    x2k = importdata ('DATA/data_AllNetworkStructures_013007.txt');
    x1l = importdata('DATA/data_112607.txt');
    x2l = importdata('DATA/data_AllNetworkStructures_112607.txt');
    x1m = importdata('DATA/data_041708.txt');
    x2m = importdata('DATA/data_AllNetworkStructures_041708.txt');
    x1n = importdata('DATA/data_042908.txt');
    x2n = importdata('DATA/data_AllNetworkStructures_042908.txt');
    x1o = importdata('DATA/data_043008.txt');
    x2o = importdata('DATA/data_AllNetworkStructures_043008.txt');
    x1p = importdata('DATA/data_043108.txt');
    x2p = importdata('DATA/data_AllNetworkStructures_043108.txt');
    x1q = importdata('DATA/data_043208.txt');
    x2q = importdata('DATA/data_AllNetworkStructures_043208.txt');
    x1r = importdata('DATA/data_052008.txt');
    x2r = importdata('DATA/data_AllNetworkStructures_052008.txt');
    x1s = importdata('DATA/data_052108.txt');
    x2s = importdata('DATA/data_AllNetworkStructures_052108.txt');
    x1t = importdata ('DATA/data_033009.txt');
    x2t = importdata ('DATA/data_AllNetworkStructures_033009.txt');
    x1u = importdata ('DATA/data_033109.txt');
    x2u = importdata ('DATA/data_AllNetworkStructures_033109.txt');
    x1v = importdata ('DATA/data_051109.txt');
    x2v = importdata ('DATA/data_AllNetworkStructures_051109.txt');

    x1ALL = [x10; x1a; x1b; x1c; x1d; x1e; x1g; x1h; x1i; x1j; x1k; x1l; x1m; x1n; x1o; x1p; x1q; x1r; x1s; x1t; x1u; x1v];
    x2ALL = [x20; x2a; x2b; x2c; x2d; x2e; x2g; x2h; x2i; x2j; x2k; x2l; x2m; x2n; x2o; x2p; x2q; x2r; x2s; x2t; x2u; x2v];
    save XALL x1ALL x2ALL;
end;

load XALL;
size(x1ALL);
size(x2ALL);
TotalObs = size(x1ALL,1);
TotalMkts = TotalObs / 4;
nMonteDraws = 400;
SegmentSize = 5540;

randn('state',0);
SDMultipleCost = .25^.5;
SDMultiplePop = .05^.5;
SD1 = 9;
SD2 = 300;
SD_HospMeasurementErrors = SDMultipleCost*SD1;
SD_NBarMeasurementErrors = SDMultiplePop*SD2;

HospMeasurementErrorsAll = normrnd(0,SD_HospMeasurementErrors, TotalObs/2, 1);
NBarMeasurementErrorsAll = normrnd(0,SD_NBarMeasurementErrors, TotalObs/4, 1);

NumHosp = 2;
NumHMO = 2;
gamma=.075;

Segment=1;
x1 = x1ALL;
x2 = x2ALL;
clear x1ALL x2ALL;

%Flags
RobustIneq=0;
IgnoreStructErrors = 0;         %Mis-specifies error by not subtracting it out
doVar = 0;                      %Calculate 95 conf interval
doTestStatistic=0;
DimBeta=1;
nConfInterval_Iterations=200;
ConfInterval = .95;

%INSTRUMENTS        DEFAULT IS FULL INSTRUMENTS
useConstINST = 0;
usemeaserrorsINST = 0;
useNoCostINST = 0;

%MEASUREMENT ERROR
usemeaserrors = 0;              %Add Measurement Noise to Moments
nMeasErrorIter = 0;
SDMultipleCost = .25^.5;
SDMultiplePop = .05^.5;
SDMultipleCost = 0;
SDMultiplePop = 0;
SD1 = 9;
SD2 = 300;
SD_HospMeasurementErrors = SDMultipleCost*SD1;
SD_NBarMeasurementErrors = SDMultiplePop*SD2;

%Probability Parameters
doProbabilities=0;
doProbBootstrap=0;
doProbVar=0;
doProbMeasErrors=1;    
doConfInterval=0;           %Conf Intervals for Prob/Meas Error case (redrawing errors, keeping data same)
nBetaProbIterations=1;
nMeasErrorIterProb=1;
nSetsOfErrors=10;
Min=16.4;
Max=17.0;
Grid=20;

%Measurement Errors
SD1 = 9;
SD2 = 300;
SD_HospMeasurementErrors = SDMultipleCost*SD1;
SD_NBarMeasurementErrors = SDMultiplePop*SD2;

if ~RobustIneq
    NumMoments = 5;             %A1: HMOs can reverse choice with largest hosp, A2: HMOs can reverse choice with smallest hosp, B:  Hosp can reverse choice with HMOs that acc
else
    NumMoments = 2;
end;

lb_temp=-100;
ub_temp=10000;
lb=lb_temp*ones(DimBeta,1);
ub=ub_temp*ones(DimBeta,1);

options = optimset('MaxFunEvals',10000, 'MaxIter',10000);

%Starting Points for Search Routines
Beta_Start = [-100, 0, 10, 100, 1000, 10000];

NumHosp = 2;
NumHMO = 2;
gamma=.075;
% Import Data %
x = [x1];
nObs = size(x,[1]);
x = [zeros(nObs,1) x];
% Prepare Data Matrices %
ID_RunID = x(:,1);
ID_Iter = x(:,2);
ID_HMO = x(:,3);
ID_Hosp = x(:,4);
NBar = x(:,5);
ISCONT = x(:,6);
HospCost_j = x(:,7);
HospCost_j2 = x(:,8);
ISCAP_j = x(:,9);
ISCAP_j2 = x(:,10);
HospCap_j = x(:,11);
HospCap_j2 = x(:,12);
HospChar = x(:,13);
Tjk = x(:,14);
Tjk2 = x(:,15);
Tj2k = x(:,16);
Njk = x(:,17);
Njk2 = x(:,18);
Nj2k = x(:,19);
HospRealProf = x(:,20);
HospRealCosts = x(:,21);
HospNewRealProf = x(:,22);
HospNewRealCosts = x(:,23);
HMONewPrem = x(:,24);
N2jk = x(:,25);
N2jk2 = x(:,26);
HMOChar = x(:,27);
HMOPrem = x(:,28);
HMOCost = x(:,29);
HMORealProf = x(:,30);
SigmaM_k = x(:,31);
SigmaM_k2 = x(:,32);
HMONewProf = x(:,33);
Sigma2M_k = x(:,34);
Sigma2M_k2 = x(:,35);
N2jk_temp = x(:,36);
N2j2k = x(:,37);
MktStd_k = x(:,38);
MktStd_k2 = x(:,39);

clear x x1;    

HospCap_Mkt = HospCap_j+ HospCap_j2;
HospCost_j_Wtd = HospCost_j.*(HospCap_j./HospCap_Mkt);
HospCost_j2_Wtd = HospCost_j2.*(HospCap_j2./HospCap_Mkt);
HospCost_MktAvg_Wtd = (HospCost_j_Wtd +HospCost_j2_Wtd)/2;
HospCost_MktAvg = (HospCost_j + HospCost_j2)/2;
HospCost_DiffAvg = HospCost_j - HospCost_MktAvg;
HospCost_DiffAvgWts = HospCost_j - HospCost_MktAvg_Wtd;
HospCost_DiffAvgWts_j2 = HospCost_j2 - HospCost_MktAvg_Wtd;

%PopPerBed = NBar./(HospCap_j + HospCap_j2);

SigmaH_jk = zeros(nObs,1);
SigmaH_j2k = zeros(nObs,1);
SigmaH_jk2 = zeros(nObs,1);
Sigma2H_jk = zeros(nObs,1);
Sigma2H_j2k = zeros(nObs,1);
Sigma2H_jk2 = zeros(nObs,1);
for idx=1:nObs
    SigmaH_jk(idx,1) = Njk(idx,1)/(SigmaM_k(idx,1)*NBar(idx,1)*gamma);
    SigmaH_j2k(idx,1) = Nj2k(idx,1)/(SigmaM_k(idx,1)*NBar(idx,1)*gamma);
    SigmaH_jk2(idx,1) = Njk2(idx,1)/(SigmaM_k2(idx,1)*NBar(idx,1)*gamma);
    Sigma2H_jk(idx,1) = N2jk(idx,1)/(Sigma2M_k(idx,1)*NBar(idx,1)*gamma);
    Sigma2H_j2k(idx,1) = N2j2k(idx,1)/(Sigma2M_k(idx,1)*NBar(idx,1)*gamma);
    Sigma2H_jk2(idx,1) = N2jk2(idx,1)/(Sigma2M_k2(idx,1)*NBar(idx,1)*gamma);
end;
DeltaNjk = Njk-N2jk;
DeltaNj2k = Nj2k-N2j2k;
DeltaNjk2 = Njk2-N2jk2;

nMkts = nObs / 4;
nContracts = 0;
for i=1:size(ISCONT,[1])
    if ISCONT(i,1)==1
        nContracts = nContracts+1;
    end;
end;
HClassTemp = [1;1;2;2;];
HospClass = kron(nMkts,HClassTemp);

ISCONTj2 = ISCONT;
ISCONTk2 = ISCONT;
for idx=1:nObs
    if (mod(idx-1,4)<2)
        ISCONTj2(idx,1) = ISCONT(idx+2,1);
    else
        ISCONTj2(idx,1) = ISCONT(idx-2,1);
    end;
    if (mod(idx,2)==1)
        ISCONTk2(idx,1) = ISCONT(idx+1,1);
    else
        ISCONTk2(idx,1) = ISCONT(idx-1,1);
    end;    
end;
DeltaISCONT = zeros(nObs,1);
for i=1:nObs
    DeltaISCONT(i,1) = (-1)^(ISCONT(i,1)+1);
    DeltaISCONTj2(i,1) = (-1)^(ISCONTj2(i,1)+1);
    DeltaISCONTk2(i,1) = (-1)^(ISCONTk2(i,1)+1);
end;

HospCostPerPatient = zeros(nObs,1);
HospCostPerPatient2 = zeros(nObs,1);
for i=1:nObs
    if (Njk(i,1) + Njk2(i,1))>0
        HospCostPerPatient(i,1) = HospRealCosts(i,1) / (Njk(i,1) + Njk2(i,1));
    end;
    if (N2jk(i,1) + N2jk2(i,1))>0    
        HospCostPerPatient2(i,1) = HospNewRealCosts(i,1) / (N2jk(i,1) + N2jk2(i,1));
    end;
end;


%==========================================================================
%INSTRUMENTS =========================
%==========================================================================
NormInst = [ones(nObs,1) Njk NBar HospCap_j HospCost_j ISCAP_j MktStd_k (HospCap_j+HospCap_j2) HMOChar HospChar];
NjkInteraction = [Njk.*Njk Njk.*NBar Njk.*HospCap_j Njk.*HospCost_j Njk.*ISCAP_j Njk.*MktStd_k Njk.*(HospCap_j+HospCap_j2) Njk.*HMOChar Njk.*HospChar];
Njk2Interaction = [Njk.*Njk.*Njk Njk.*Njk.*NBar Njk.*Njk.*HospCap_j Njk.*Njk.*HospCost_j Njk.*Njk.*ISCAP_j Njk.*Njk.*MktStd_k Njk.*Njk.*(HospCap_j+HospCap_j2) Njk.*Njk.*HMOChar Njk.*Njk.*HospChar];
MktInteraction = [MktStd_k.*NBar MktStd_k.*HospCap_j MktStd_k.*HospCost_j MktStd_k.*ISCAP_j MktStd_k.*MktStd_k MktStd_k.*(HospCap_j+HospCap_j2) MktStd_k.*HMOChar MktStd_k.*HospChar];
ISCAPInteraction = [ISCAP_j.*NBar ISCAP_j.*HospCap_j ISCAP_j.*HospCost_j ISCAP_j.*(HospCap_j+HospCap_j2) ISCAP_j.*HMOChar ISCAP_j.*HospChar];
CostInteraction = [HospCost_j.*NBar HospCost_j.*HospCap_j HospCost_j.*HospCost_j HospCost_j.*(HospCap_j+HospCap_j2) HospCost_j.*HMOChar HospCost_j.*HospChar];
SumHospCapInteraction = [HospCap_j.*(HospCap_j+HospCap_j2) HMOChar.*(HospCap_j+HospCap_j2) HospChar.*(HospCap_j+HospCap_j2) (HospCap_j+HospCap_j2).*(HospCap_j+HospCap_j2) ];
CharInteraction = [HospChar.*HospChar HospChar.*HMOChar HMOChar.*HMOChar];
HospCostsWtd = [HospCost_j_Wtd HospCost_j2_Wtd HospCost_MktAvg HospCost_MktAvg_Wtd];

INST = [NormInst NjkInteraction Njk2Interaction MktInteraction ISCAPInteraction CostInteraction SumHospCapInteraction CharInteraction HospCostsWtd];
INST_NoContract = [ones(nObs,1) NBar HospCap_j HospCost_j ISCAP_j MktStd_k (HospCap_j+HospCap_j2) HMOChar HospChar MktInteraction ISCAPInteraction CostInteraction SumHospCapInteraction CharInteraction HospCostsWtd];

ALLINST_NoNjk = INST;
ALLINST_NoNjk(:,2) = [];

NBarMkt = zeros(nMkts,1);
SumCapMkt = zeros(nMkts,1);
SumCostMkt = zeros(nMkts,1);
SumISCAPMkt = zeros(nMkts,1);
for iMkt=1:nMkts
    idx = (iMkt-1)*4+1;
    NBarMkt(iMkt,1) = NBar(idx,1);
    SumCapMkt(iMkt,1) = (HospCap_j(idx,1)+HospCap_j2(idx,1));
    SumCostMkt(iMkt,1) = (HospCost_j(idx,1)+HospCost_j2(idx,1));
    SumISCAPMkt(iMkt,1) = (ISCAP_j(idx,1)+ISCAP_j2(idx,1));    
end;

MktINST = [ones(nMkts,1) NBarMkt SumCapMkt SumCostMkt SumISCAPMkt SumCapMkt./NBarMkt SumCostMkt./SumCapMkt];
MktINSTmeas = MktINST;
if ~(SDMultiplePop + SDMultipleCost == 0)
    if (SDMultiplePop>0) && (SDMultipleCost>0)
        %pop and cost measurement error
        MktINSTmeas = [ones(nMkts,1) SumCapMkt SumISCAPMkt];
    else
        %Cost only measurement error
        MktINSTmeas = [ones(nMkts,1) NBarMkt SumCapMkt SumISCAPMkt SumCapMkt./NBarMkt];
    end;
end; 
if useConstINST
    INST = [ones(nObs,1) ];
    INST_NoContract = [ones(nObs,1)];
    MktINST = [ones(nMkts,1)];
    MktINSTmeas = [ones(nMkts,1)];
end;

nInst = size(INST,2);
nInstNoContract = size(INST_NoContract,2);




%==========================================================================
%ERROR CALCULATIONS========================================================
%==========================================================================
                        useTotalTransfers = 1;
                        %Calculate distribution of errors using total transfers (need to account
                        %for Njk/N2jk issue
                        if useTotalTransfers
                            NNjk=Njk+N2jk;
                            NNj2k=NNjk;
                            NNjk2=NNjk;
                            NN2jk=NNjk;
                            NN2j2k=NNjk;
                            NN2jk2=NNjk;
                            for idx=1:nObs
                                if (mod(idx-1,4)<2)
                                    NNj2k(idx,1) = NNjk(idx+2,1);
                                else
                                    NNj2k(idx,1) = NNjk(idx-2,1);
                                end;
                                if (mod(idx,2)==1)
                                    NNjk2(idx,1) = NNjk(idx+1,1);
                                else
                                    NNjk2(idx,1) = NNjk(idx-1,1);
                                end;    
                            end;

                            X = [NNjk ALLINST_NoNjk];
                            Beta = (inv(X'*X)*(X'*(NNjk.*Tjk)));
                            Beta_0 = Beta(1);
                            TrueBeta = Beta;
                            TrueBeta0 = Beta_0;
                            Tjk_NEW = (NNjk.*Tjk - X*TrueBeta + NNjk*TrueBeta0)./NNjk;

                            Tj2k_NEW = Tjk_NEW;
                            Tjk2_NEW = Tjk_NEW; 

                            Tjk_Residual = Tjk - Tjk_NEW;
                            Tj2k_Residual = Tjk_Residual;
                            Tjk2_Residual = Tjk_Residual;    

                            for idx=1:nObs
                                if (mod(idx-1,4)<2)
                                    Tj2k_NEW(idx,1) = Tjk_NEW(idx+2,1);
                                    Tj2k_Residual(idx,1) = Tjk_Residual(idx+2,1);
                                else
                                    Tj2k_NEW(idx,1) = Tjk_NEW(idx-2,1);
                                    Tj2k_Residual(idx,1) = Tjk_Residual(idx-2,1);
                                end;
                                if (mod(idx,2)==1)
                                    Tjk2_NEW(idx,1) = Tjk_NEW(idx+1,1);
                                    Tjk2_Residual(idx,1) = Tjk_Residual(idx+1,1);
                                else
                                    Tjk2_NEW(idx,1) = Tjk_NEW(idx-1,1);
                                    Tjk2_Residual(idx,1) = Tjk_Residual(idx-1,1);
                                end;    
                            end;  

                            Ejk = Njk.* Tjk_NEW - Njk * Beta_0;    
                            Ej2k = Nj2k.*Tj2k_NEW - Nj2k * Beta_0;    
                            Ejk2 = Njk2.*Tjk2_NEW - Njk2 * Beta_0;
                            E2jk = N2jk.*Tjk_NEW - N2jk * Beta_0;
                            E2j2k = N2j2k.*Tj2k_NEW - N2j2k * Beta_0;    
                            E2jk2 = N2jk2.*Tjk2_NEW - N2jk2 * Beta_0;   

                            Delta_Ejk = Ejk - E2jk;
                            Delta_Ej2k = Ej2k - E2j2k;
                            Delta_Ejk2 = Ejk2 - E2jk2;
                        else
                            cons = ones(nObs,1);
                            Beta_0 = (inv(cons'*cons)*(cons'*(Tjk)));        
                            Ejk = Tjk - Beta_0;  
                            Ej2k = Tj2k - Beta_0;
                            Ejk2 = Tjk2 - Beta_0;
                        end;

                        Err2=zeros(nMkts,4);                % Residuals in nMkts x 4 form. Reshaped Err.
                        for i=0:nMkts-1
                            for j=1:4
                                Err2(i+1,j) = Ejk(i*4+j);
                            end;
                        end;
                        ErrMean = mean(Err2)';
                        VarErr=zeros(4,4);
                        for n=1:nMkts
                            SSample = Err2(n,:)';
                            VarErr = VarErr + ((SSample-ErrMean)*(SSample-ErrMean)');
                        end;
                        VarErr = VarErr / (nMkts);
                        CholVarErr = chol(VarErr);        

                        %display(Beta_0);
                        Beta_Start = [Beta_Start Beta_0];
                        NumBetaStart = size(Beta_Start,2);   

                    RealNu2s = Err2;
clear X Ejk Ej2k Ejk2;
save Estimates_TrueBetaNu2s_072709.mat CholVarErr RealNu2s TrueBeta0 ErrMean TrueBeta VarErr;

%==========================================================================
%PRELIMINARY MOMENT CONSTRUCTIONS====================================================
%==========================================================================
    DeltaHospCosts = HospRealCosts - HospNewRealCosts;
    DeltaHMOPrem = ( SigmaM_k.*(HMOPrem-HMOCost) - Sigma2M_k.*(HMONewPrem-HMOCost)) .*NBar;

    AbsDeltaNjkInv = abs(ones(nObs,1) ./ DeltaNjk);

    A_Hosp = DeltaNjk + DeltaNjk2;
    A_HMO = (-1) * (DeltaNjk + DeltaNj2k);
    C_Hosp = DeltaHospCosts - DeltaNjk.*Tjk_Residual - DeltaNjk2.*Tjk2_Residual - Delta_Ejk2;
    C_HMO = (-1)*(DeltaHMOPrem - DeltaNjk.*Tjk_Residual - DeltaNj2k.*Tj2k_Residual - Delta_Ej2k);

    if (~IgnoreStructErrors)
        C_Hosp = C_Hosp - Delta_Ejk;
        C_HMO = C_HMO - (-1)*(Delta_Ejk);          %Error term is "added" to C_HMO, which is subtracted from profits
    end;

    A_HMO_Orig = A_HMO;
    A_Hosp_Orig = A_Hosp;
    C_HMO_Orig = C_HMO;
    C_Hosp_Orig = C_Hosp;

    DeltaHospTransReal = DeltaNjk.*Tjk+DeltaNjk2.*Tjk2;
    DeltaHMOTransReal =  DeltaNjk.*Tjk+DeltaNj2k.*Tj2k;
    DeltaHospTransPredicted = (DeltaNjk+DeltaNjk2)*Beta_0;
    DeltaHMOTransPredicted =  (DeltaNjk+DeltaNj2k)*Beta_0;

%==========================================================================
%BEGIN CONSTRUCTION OF MOMENTS====================================================
%==========================================================================
if usemeaserrors
    MeasErrorBetas = zeros(nMeasErrorIter,2);
    MeasErrorSinglePoints=0;
else
    nMeasErrorIter = 1;
end;

for idxMeasError=1:nMeasErrorIter
    if usemeaserrors
        HospMeasurementErrors = normrnd(0,SD_HospMeasurementErrors, nObs/2, 1);
        Etaj = kron(HospMeasurementErrors,ones(2,1));
        NBarMeasurementErrors = normrnd(0,SD_NBarMeasurementErrors, nObs/4, 1);
        Omega = kron(NBarMeasurementErrors,ones(4,1));    

        Psijk_A = gamma * (SigmaM_k.*SigmaH_jk);
        Psijk_B = gamma * (Sigma2M_k.*Sigma2H_jk);
        Psijk = Psijk_A-Psijk_B;
        Psij2k_A = gamma * (SigmaM_k.*SigmaH_j2k);
        Psij2k_B = gamma * (Sigma2M_k.*Sigma2H_j2k);
        Psij2k = Psij2k_A-Psij2k_B;        
        Psijk2_A = gamma * (SigmaM_k2.*SigmaH_jk2);
        Psijk2_B = gamma * (Sigma2M_k2.*Sigma2H_jk2);
        Psijk2 = Psijk2_A-Psijk2_B;        

        PsiHosp_A = Psijk_A + Psijk2_A;
        PsiHosp_B = Psijk_B + Psijk2_B;

        DeltaNjkErr = DeltaNjk + Omega.*Psijk;
        AbsDeltaNjkErrInv = abs(ones(nObs,1) ./ DeltaNjkErr);
        DeltaNjk2Err = DeltaNjk2 + Omega.*Psijk2;
        DeltaNj2kErr = DeltaNj2k + Omega.*Psij2k;        

        A_Hosp = DeltaNjk + DeltaNjk2;
        A_HMO = (-1) * (DeltaNjk + DeltaNj2k);
        C_Hosp = (Njk+Njk2+Omega.*PsiHosp_A).*(HospCostPerPatient+Etaj) - (N2jk+N2jk2+Omega.*PsiHosp_B).*(HospCostPerPatient2+Etaj) - DeltaNjk.*Tjk_Residual - DeltaNjk2.*Tjk2_Residual - Delta_Ejk2;
        C_HMO = (-1)*(   ( SigmaM_k.*(HMOPrem-HMOCost)-Sigma2M_k.*(HMONewPrem-HMOCost) ).*(NBar+Omega) - DeltaNjk.*Tjk_Residual - DeltaNj2k.*Tj2k_Residual- Delta_Ej2k );

        if (~IgnoreStructErrors)
            C_Hosp = C_Hosp - Delta_Ejk;
            C_HMO = C_HMO - (-1)*(Delta_Ejk);          %Error term is "added" to C_HMO, which is subtracted from profits
        end;

    end;

    %Create Moment Matrices (normal case)=======================================================================
    if ~RobustIneq
        Z_AB = zeros(nInst,DimBeta,nMkts,NumMoments);
        W_AB = zeros(nInst,1,nMkts,NumMoments);
        for MomentID=1:NumMoments
            % Moments for HMO: If HMO Accepts, should be better off. Moment
            % A1 for larger hospital; Moment A2 for smaller hospital.
            if (MomentID<3)
                for j=1:nMkts
                    for k=1:NumHMO
                        idx = (j-1)*(NumHosp*NumHMO) + (MomentID-1)*(NumHosp) + k;
                        hTemp = INST(idx,:)';
                        hTempNoContract = INST_NoContract(idx,:)';
                        Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID) + kron( A_HMO(idx,:) , hTemp);
                        W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID) + kron( C_HMO(idx,:) , hTemp);
                    end;
                    Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID)/ 4;
                    W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID)/ 4;
                end;
            else
                % Moments for Hospitals: Hospitals that do contract should be
                % better off than if they did not contract
                if (MomentID==3)
                    for j=1:nMkts
                        Count = 0;
                        for k=1:(NumHosp*NumHMO)
                            idx = (j-1)*(NumHosp*NumHMO) + k;
                            hTemp = INST(idx,:)';
                            if ISCONT(idx,1)==1
                                Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID) + kron( A_Hosp(idx,:) , hTemp);
                                W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID) + kron( C_Hosp(idx,:) , hTemp);
                                Count = Count + 1;
                            end;
                        end;
                        if Count==0
                            Count = 1;
                        end;
                        Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID)/ Count;
                        W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID)/ Count;                    
                    end;    
                end;%ifMomentID==5
            end;%ifMomentID<5     
            if (MomentID > 3)
                if (MomentID==4)
                    for j=1:nMkts
                        for k=1:(NumHosp*NumHMO)
                           idx = (j-1)*(NumHosp*NumHMO) + k;
                           hTemp = INST(idx,:)';
                           Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID) + kron( ISCONT(idx,1)*A_Hosp(idx,:) + (1-ISCONT(idx,1))*A_HMO(idx,:) , hTemp) / (NumHosp*NumHMO);
                           W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID) + kron( ISCONT(idx,1)*C_Hosp(idx,:) + (1-ISCONT(idx,1))*C_HMO(idx,:) , hTemp) / (NumHosp*NumHMO);
                        end;
                    end;
                end;
                if MomentID==5
                    if NumMoments >= 2    
                        % Sum of Proift changes if contract should be positive if they
                        % do contract
                        for j=1:nMkts
                            Count = 0;
                            for k=1:(NumHosp*NumHMO)
                                idx = (j-1)*(NumHosp*NumHMO) + k;
                                hTemp = INST(idx,:)';
                                Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID) + kron( (ISCONT(idx,1)).*(A_Hosp(idx,:)+A_HMO(idx,:)) , hTemp);
                                W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID) + kron( (ISCONT(idx,1)).*(C_Hosp(idx,:)+C_HMO(idx,:)) , hTemp);
                                if ISCONT(idx,1)==1 
                                    Count = Count + 1;
                                end;
                            end;
                            if Count==0
                                Count = 1;
                            end;
                            Z_AB(:,:,j,MomentID) = Z_AB(:,:,j,MomentID)/ Count;
                            W_AB(:,:,j,MomentID) = W_AB(:,:,j,MomentID)/ Count;                       
                        end;
                    end;   
                end;
            end;                                 
        end;%forMomentID=1:NumMoments
    else    %==========================================COMPUTE ROBUST INEQUALITIES: ===========================
        Z_AB = zeros(nInst,DimBeta,nMkts,NumMoments);
        W_AB = zeros(nInst,1,nMkts,NumMoments);

        % Adding change in profits for hospital if do contract and 
        % change in profits for HMO if don't contract
        for j=1:nMkts
            for k=1:(NumHosp*NumHMO)
               idx = (j-1)*(NumHosp*NumHMO) + k;
               hTemp = INST(idx,:)';
               Z_AB(:,:,j,1) = Z_AB(:,:,j,1) + kron( ISCONT(idx,1)*A_Hosp(idx,:) + (1-ISCONT(idx,1))*A_HMO(idx,:) , hTemp) / (NumHosp*NumHMO);
               W_AB(:,:,j,1) = W_AB(:,:,j,1) + kron( ISCONT(idx,1)*C_Hosp(idx,:) + (1-ISCONT(idx,1))*C_HMO(idx,:) , hTemp) / (NumHosp*NumHMO);
            end;
        end;

        if NumMoments >= 2    
            % Sum of Proift changes if contract should be positive if they
            % do contract
            for j=1:nMkts
                Count = 0;
                for k=1:(NumHosp*NumHMO)
                    idx = (j-1)*(NumHosp*NumHMO) + k;
                    hTemp = INST(idx,:)';
                    Z_AB(:,:,j,2) = Z_AB(:,:,j,2) + kron( (ISCONT(idx,1)).*(A_Hosp(idx,:)+A_HMO(idx,:)) , hTemp);
                    W_AB(:,:,j,2) = W_AB(:,:,j,2) + kron( (ISCONT(idx,1)).*(C_Hosp(idx,:)+C_HMO(idx,:)) , hTemp);
                    if ISCONT(idx,1)==1 
                        Count = Count + 1;
                    end;
                end;
                if Count==0
                    Count = 1;
                end;
                Z_AB(:,:,j,2) = Z_AB(:,:,j,2)/ Count;
                W_AB(:,:,j,2) = W_AB(:,:,j,2)/ Count;                       
            end;
        end;
    end;    %if RobustIneq
    %==========================================================================
    ZJtempAB = zeros(nInst,DimBeta,NumMoments);
    WJtempAB = zeros(nInst,1,NumMoments);
    for j=1:nMkts
        for MomentID=1:NumMoments
            ZJtempAB(:,:,MomentID) = ZJtempAB(:,:,MomentID) + Z_AB(:,:,j,MomentID)/nMkts;
            WJtempAB(:,1,MomentID) = WJtempAB(:,1,MomentID) + W_AB(:,1,j,MomentID)/nMkts;
        end;
    end;
    ZJ =[];
    WJ =[];
    for n=1:NumMoments
        ZJ = [ZJ;ZJtempAB(:,:,n)];
        WJ = [WJ;WJtempAB(:,:,n)];
    end;

    %Testing====================================================================================================    
    SinglePoint = 0;
    Beta_S1=Beta_0;
    [BetaMin BetaMax ErrorFlag] = f_fmincon_092605(ZJ,WJ,lb,ub,options,Beta_Start);
    if (ErrorFlag == 0)
        TEMP = [Beta_0 BetaMin BetaMax];
    else
        [BetaPoint fval ErrorFlag2] = f_fminconsingle(ZJ,WJ,lb,ub,options,Beta_Start);
        TEMP = [Beta_0 BetaPoint BetaPoint];
        display(fval);
        TJTestStatistic = (nMkts^.5) * (fval);
        SinglePoint=1;
    end;
    display(idxMeasError);
    BETA_ID = TEMP;
    display(BETA_ID);

    if usemeaserrors
        MeasErrorBetas(idxMeasError,1) = TEMP(1,2);
        MeasErrorBetas(idxMeasError,2) = TEMP(1,3);    
        if SinglePoint
            MeasErrorSinglePoints = MeasErrorSinglePoints+1;
        end;
    end;
end; %END MEASUREMENT ERROR ITERATIONS
if usemeaserrors&&(nMeasErrorIter>1)        
    MeasErrorConfInterval = zeros(3,2);
    BMinMeas = sortrows(MeasErrorBetas(:,1));
    BMaxMeas = sortrows(MeasErrorBetas(:,2));
    MeasErrorConfInterval(2,1)=BMinMeas( .025*nMeasErrorIter+1,1 );
    MeasErrorConfInterval(3,1)=BMinMeas( .975*nMeasErrorIter,1 );
    MeasErrorConfInterval(2,2)=BMaxMeas( .025*nMeasErrorIter+1,1 );
    MeasErrorConfInterval(3,2)=BMaxMeas( .975*nMeasErrorIter,1 );        
    MeasErrorConfInterval(1,1)=mean(BMinMeas);
    MeasErrorConfInterval(1,2)=mean(BMaxMeas);
    display(MeasErrorConfInterval);
end;

%==========================================================================
%Get Variance on Bounds====================================================
%==========================================================================
if (doVar)
    SizeMoments = size(ZJ,1);
    SJ = [reshape(ZJ,SizeMoments*DimBeta,1); WJ];
    VarSJ = zeros(SizeMoments*(DimBeta+1),SizeMoments*(DimBeta+1));
    SSample = zeros(SizeMoments*(DimBeta+1),1);
    ZJ_Simul_All = zeros(SizeMoments,nConfInterval_Iterations);
    WJ_Simul_All = zeros(SizeMoments,nConfInterval_Iterations);
    for j=1:nMkts
        ZSample =[];
        WSample =[];
        for n=1:NumMoments
            ZSample = [ZSample;Z_AB(:,:,j,n)];
            WSample = [WSample;W_AB(:,:,j,n)];
            if (MomentID<3)&&(~RobustIneq)            
                ZSample = [ZSample;Z_C(:,:,j,n)];
                WSample = [WSample;W_C(:,:,j,n)];                
            end;
        end;
        SSample = [reshape(ZSample,SizeMoments*DimBeta,1); WSample];
        VarSJ = VarSJ + ((SSample-SJ)*(SSample-SJ)')/nMkts;
    end;
    VarSJ = (.99999)*VarSJ + (.00001)*eye(size(VarSJ,1),size(VarSJ,1));
    BetaSimulationMin = zeros(size(Beta_0,1),nConfInterval_Iterations);
    BetaSimulationMax = zeros(size(Beta_0,1),nConfInterval_Iterations);
    temp = chol(VarSJ);

    NumSinglePoints = 0;

    for n=1:nConfInterval_Iterations
        epsilons = normrnd(0,1, SizeMoments*(DimBeta+1), 1);
        UJ_simul = temp' * epsilons / ((nMkts)^.5);
        SJ_Simul = SJ + UJ_simul;
        ZJWJ_Simul = reshape(SJ_Simul,SizeMoments,DimBeta+1);
        ZJ_Simul = ZJWJ_Simul(:,1:DimBeta);
        WJ_Simul = ZJWJ_Simul(:,DimBeta+1);

        ZJ_Simul_All(:,n) = UJ_simul(1:SizeMoments,1) * ((nMkts)^.5);
        WJ_Simul_All(:,n) = UJ_simul(SizeMoments+1:SizeMoments*2,1) * ((nMkts)^.5);

        [BetaSimulationMinTemp,BetaSimulationMaxTemp, ErrorFlag] = f_fmincon_092605(ZJ_Simul,WJ_Simul,lb,ub,options,Beta_Start);
        if (ErrorFlag == 0)
            BetaSimulationMin(:,n) = BetaSimulationMinTemp;
            BetaSimulationMax(:,n) = BetaSimulationMaxTemp;
        else
            [BetaPointTemp fval ErrorFlag] = f_fminconsingle(ZJ_Simul,WJ_Simul,lb,ub,options,Beta_Start);
            BetaSimulationMin(:,n) = BetaPointTemp;
            BetaSimulationMax(:,n) = BetaPointTemp;
            NumSinglePoints=NumSinglePoints+1;
        end;
        display(n);
    end;
    BetaConfInterval = zeros(3,2);

    BMin = sortrows(BetaSimulationMin(1,:)');
    BMax = sortrows(BetaSimulationMax(1,:)');
    BetaConfInterval(1,1)=BMin( .025*nConfInterval_Iterations+1,1 );
    BetaConfInterval(2,1)=BMin( .025*nConfInterval_Iterations+1,1 );
    BetaConfInterval(3,1)=BMin( .975*nConfInterval_Iterations,1 );
    BetaConfInterval(1,2)=BMax( .975*nConfInterval_Iterations,1 );        
    BetaConfInterval(2,2)=BMax( .025*nConfInterval_Iterations+1,1 );
    BetaConfInterval(3,2)=BMax( .975*nConfInterval_Iterations,1 );        

    BMinMax = [BetaSimulationMin(1,:)' BetaSimulationMax(1,:)'];
    display(BetaConfInterval);
    display(NumSinglePoints);


    if doTestStatistic
        BetaConfInterval2 = zeros(1,2);
        BetaConfInterval2(1,1) = BMin( .0125*nConfInterval_Iterations+1,1 ); 
        BetaConfInterval2(1,2) = BMax( .9875*nConfInterval_Iterations,1 );
        [zJTestStatistic BetaActual zJ0TestStatistic] = f_zJTestStatistic(BetaConfInterval2, ZJ_Simul_All, WJ_Simul_All, BetaPoint);
        TestStatistic = TJTestStatistic/zJTestStatistic;
        TestStatistic0 = TJTestStatistic/zJ0TestStatistic;    
        display(TestStatistic);
        display(TestStatistic0);    
    end;
end; %if doVar

%==========================================================================
%==========================================================================
%==========================================================================
% PROBABILITIES
%==========================================================================
%==========================================================================
%==========================================================================
if doProbabilities
    ISCONT_All = zeros(nObs,16);
    Njk_All = zeros(nObs,16);
    Njk2_All = zeros(nObs,16);
    Nj2k_All = zeros(nObs,16);
    N2jk_All = zeros(nObs,16);
    N2jk2_All = zeros(nObs,16);
    N2j2k2_All = zeros(nObs,16);
    Prem_All = zeros(nObs,16);
    Prem2_All = zeros(nObs,16);
    HospCost_All = zeros(nObs,16);
    HospCost2_All = zeros(nObs,16);
    SigmaM_k_All = zeros(nObs,16);    
    SigmaM_k2_All = zeros(nObs,16);     
    SigmaM2_k_All = zeros(nObs,16);
    SigmaM2_k2_All = zeros(nObs,16);    
    HospCost_j_All =  zeros(nObs,16);    
    HospCost_j2_All =  zeros(nObs,16);        
    RealHospProf_All = zeros(nObs,16);    
    RealHospProf2_All = zeros(nObs,16);
    RealHospCost_All = zeros(nObs,16);
    RealHospCost2_All = zeros(nObs,16);
    RealHMOProf_all = zeros(nObs,16);
    RealHMOProf2_all = zeros(nObs,16);

    DeltaNjk_All = zeros(nObs,16);    
    DeltaNj2k_All = zeros(nObs,16);    
    DeltaNjk2_All = zeros(nObs,16);        
    DeltaNjkInv_All = zeros(nObs,16);        

    ISCONTj2_All = zeros(nObs,16);
    ISCONTk2_All = zeros(nObs,16);
    DeltaISCONT_All = zeros(nObs,16);
    DeltaISCONTj2_All = zeros(nObs,16);
    DeltaISCONTk2_All = zeros(nObs,16);

    SigmaH_jk_All = zeros(nObs,16);
    SigmaH_j2k_All = zeros(nObs,16);
    SigmaH_jk2_All = zeros(nObs,16);
    Sigma2H_jk_All = zeros(nObs,16);
    Sigma2H_j2k_All = zeros(nObs,16);
    Sigma2H_jk2_All = zeros(nObs,16);    

    DeltaHMOPremTemp_All = zeros(nObs,16);

    EqNetwork = x2(:,4);
    x2Idx=9; %where alternative network structure data begins
    for i = 1:16
        ISCONT_All(:,i) = x2(:,x2Idx);
        Njk_All(:,i) = x2(:,x2Idx+1);
        Njk2_All(:,i) = x2(:,x2Idx+2);
        Nj2k_All(:,i) = x2(:,x2Idx+3);        
        N2jk_All(:,i) = x2(:,x2Idx+4);
        N2jk2_All(:,i) = x2(:,x2Idx+5);
        N2j2k_All(:,i) = x2(:,x2Idx+6);                
        Prem_All(:,i) = x2(:,x2Idx+7);
        Prem2_All(:,i) = x2(:,x2Idx+8);
        SigmaM_k_All(:,i) = x2(:,x2Idx+9);
        SigmaM_k2_All(:,i) = x2(:,x2Idx+10);
        SigmaM2_k_All(:,i) = x2(:,x2Idx+11);
        SigmaM2_k2_All(:,i) = x2(:,x2Idx+12);
        HospCost_j_All(:,i) = x2(:,x2Idx+13);
        HospCost_j2_All(:,i) = x2(:,x2Idx+14);        
        RealHospProf_All(:,i) = x2(:,x2Idx+15);
        RealHospProf2_All(:,i) = x2(:,x2Idx+16);
        RealHospCost_All(:,i) = x2(:,x2Idx+17);
        RealHospCost2_All(:,i) = x2(:,x2Idx+18);
        RealHMOProf_All(:,i) = x2(:,x2Idx+19);
        RealHMOProf_All(:,i) = x2(:,x2Idx+20);        
        x2Idx = x2Idx + 22;           

        load Estimates_TrueBetaNu2s_072709.mat            

        DeltaNjk_All(:,i) = Njk_All(:,i)-N2jk_All(:,i);
        DeltaNj2k_All(:,i) = Nj2k_All(:,i)-N2j2k_All(:,i);
        DeltaNjk2_All(:,i) = Njk2_All(:,i)-N2jk2_All(:,i);
        AbsDeltaNjkInv_All(:,i) = abs(ones(nObs,1) ./ DeltaNjk_All(:,i));

        ISCONTj2_All(:,i) = ISCONT_All(:,i);
        ISCONTk2_All(:,i) = ISCONT_All(:,i);
        for idx=1:nObs
            if (mod(idx-1,4)<2)
                ISCONTj2_All(idx,i) = ISCONT_All(idx+2,i);
            else
                ISCONTj2_All(idx,i) = ISCONT_All(idx-2,i);
            end;
            if (mod(idx,2)==1)
                ISCONTk2_All(idx,i) = ISCONT_All(idx+1,i);
            else
                ISCONTk2_All(idx,i) = ISCONT_All(idx-1,i);
            end;
        end;
        DeltaISCONT_All(:,i) = zeros(nObs,1);
        for idx=1:nObs
            DeltaISCONT_All(idx,i) = (-1)^(ISCONT_All(idx,i)+1);
            DeltaISCONTj2_All(idx,i) = (-1)^(ISCONTj2_All(idx,i)+1);
            DeltaISCONTk2_All(idx,i) = (-1)^(ISCONTk2_All(idx,i)+1);

            SigmaH_jk_All(idx,i) = Njk_All(idx,i)/(SigmaM_k_All(idx,i)*NBar(idx,1)*gamma);
            SigmaH_j2k_All(idx,i) = Nj2k_All(idx,i)/(SigmaM_k_All(idx,i)*NBar(idx,1)*gamma);
            SigmaH_jk2_All(idx,i) = Njk2_All(idx,i)/(SigmaM_k2_All(idx,i)*NBar(idx,1)*gamma);
            Sigma2H_jk_All(idx,i) = N2jk_All(idx,i)/(SigmaM2_k_All(idx,i)*NBar(idx,1)*gamma);
            Sigma2H_j2k_All(idx,i) = N2j2k_All(idx,i)/(SigmaM2_k_All(idx,i)*NBar(idx,1)*gamma);
            Sigma2H_jk2_All(idx,i) = N2jk2_All(idx,i)/(SigmaM2_k2_All(idx,i)*NBar(idx,1)*gamma);
        end;
        HospCostPerPatient_All = zeros(nObs,16);
        HospCostPerPatient2_All = zeros(nObs,16);
        for idxStruct=1:16
            for idxObs=1:nObs
                if ( (Njk_All(idxObs,idxStruct) + Njk2_All(idxObs,idxStruct))>0 )
                    HospCostPerPatient_All(idxObs,idxStruct) = RealHospCost_All(idxObs,idxStruct) / (Njk_All(idxObs,idxStruct) + Njk2_All(idxObs,idxStruct));
                end;
                if ( (N2jk_All(idxObs,idxStruct) + N2jk2_All(idxObs,idxStruct))>0 )
                    HospCostPerPatient2_All(idxObs,idxStruct) = RealHospCost2_All(idxObs,idxStruct) / (N2jk_All(idxObs,idxStruct) + N2jk2_All(idxObs,idxStruct));
                end;
            end;
        end;    
        DeltaHospCosts = HospRealCosts - HospNewRealCosts;
        DeltaHMOPremTemp_All(:,i) = ( SigmaM_k_All(:,i).*(Prem_All(:,i)-HMOCost) - SigmaM2_k_All(:,i).*(Prem2_All(:,i)-HMOCost) ).*NBar;

    end; %for i=1:16
clear x2;

    A_Hosp_All = DeltaNjk_All + DeltaNjk2_All;
    A_HMO_All = (-1) * (DeltaNjk_All + DeltaNj2k_All);
    C_Hosp_All = RealHospCost_All - RealHospCost2_All - DeltaNjk_All.*kron(Tjk_Residual,ones(1,16)) - DeltaNjk2_All.*kron(Tjk2_Residual,ones(1,16)) - kron(Delta_Ejk2,ones(1,16));
    C_HMO_All = (-1)*(DeltaHMOPremTemp_All - DeltaNjk_All.*kron(Tjk_Residual,ones(1,16)) - DeltaNj2k_All.*kron(Tj2k_Residual,ones(1,16)) - kron(Delta_Ej2k,ones(1,16)) );

    if ~IgnoreStructErrors
        C_Hosp_All = C_Hosp_All - kron(Delta_Ejk,ones(1,16));
        C_HMO_All = C_HMO_All - (-1)*(kron(Delta_Ejk,ones(1,16)));
    end;

    % 10 Cases to consider: 1-NoContracts; 2 LHosp(1),S(0); 3 L(2),S(0);
    % 4-L(0),S(1); 5-L(0),S(2); 6-L1S1 Same; 7-L1S1Diff; 8-L1S2; 9-L2S1;
    % 10-AllContract

    Struct_Indicator = zeros(nMkts,10);
    for idxMkt = 1:nMkts
        idxObs = (idxMkt-1)*4+1;
        if ( EqNetwork(idxObs,1) == 0 ) Struct_Indicator(idxMkt,1) = 1; end;
        if ( (EqNetwork(idxObs,1) == 1)||(EqNetwork(idxObs,1) == 4) ) Struct_Indicator(idxMkt,2) = 1; end;
        if ( EqNetwork(idxObs,1) == 5 ) Struct_Indicator(idxMkt,3) = 1; end;
        if ( (EqNetwork(idxObs,1) == 2)||(EqNetwork(idxObs,1) == 8) ) Struct_Indicator(idxMkt,4) = 1; end;        
        if ( EqNetwork(idxObs,1) == 10 ) Struct_Indicator(idxMkt,5) = 1; end;        
        if ( (EqNetwork(idxObs,1) == 3)||(EqNetwork(idxObs,1) == 12) ) Struct_Indicator(idxMkt,6) = 1; end;        
        if ( (EqNetwork(idxObs,1) == 6)||(EqNetwork(idxObs,1) == 9) ) Struct_Indicator(idxMkt,7) = 1; end;        
        if ( (EqNetwork(idxObs,1) == 11)||(EqNetwork(idxObs,1) == 14) ) Struct_Indicator(idxMkt,8) = 1; end;        
        if ( (EqNetwork(idxObs,1) == 7)||(EqNetwork(idxObs,1) == 13) ) Struct_Indicator(idxMkt,9) = 1; end;                
        if ( EqNetwork(idxObs,1) == 15 ) Struct_Indicator(idxMkt,10) = 1; end;                
    end;

    % Form what Structures look like in ALL NETWORK STRUCTURS
    ISCONT_AllNetworks = zeros(nObs,16);
    DeltaISCONT_AllNetworks = zeros(nObs,16);
    for idxObs = 1:nObs
        if (mod(idxObs,4)==1)
            ISCONT_AllNetworks(idxObs,5) = 1;
            ISCONT_AllNetworks(idxObs,6) = 1;
            ISCONT_AllNetworks(idxObs,7) = 1;
            ISCONT_AllNetworks(idxObs,8) = 1;
            ISCONT_AllNetworks(idxObs,13) = 1;
            ISCONT_AllNetworks(idxObs,14) = 1;
            ISCONT_AllNetworks(idxObs,15) = 1;
            ISCONT_AllNetworks(idxObs,16) = 1;            
        end;
        if (mod(idxObs,4)==2)
            ISCONT_AllNetworks(idxObs,2) = 1;
            ISCONT_AllNetworks(idxObs,4) = 1;
            ISCONT_AllNetworks(idxObs,6) = 1;
            ISCONT_AllNetworks(idxObs,8) = 1;
            ISCONT_AllNetworks(idxObs,10) = 1;
            ISCONT_AllNetworks(idxObs,12) = 1;
            ISCONT_AllNetworks(idxObs,14) = 1;
            ISCONT_AllNetworks(idxObs,16) = 1;
        end;    
        if (mod(idxObs,4)==3)
            ISCONT_AllNetworks(idxObs,9) = 1;
            ISCONT_AllNetworks(idxObs,10) = 1;
            ISCONT_AllNetworks(idxObs,11) = 1;
            ISCONT_AllNetworks(idxObs,12) = 1;
            ISCONT_AllNetworks(idxObs,13) = 1;
            ISCONT_AllNetworks(idxObs,14) = 1;
            ISCONT_AllNetworks(idxObs,15) = 1;
            ISCONT_AllNetworks(idxObs,16) = 1;            
        end;    
        if (mod(idxObs,4)==0)
            ISCONT_AllNetworks(idxObs,3) = 1;
            ISCONT_AllNetworks(idxObs,4) = 1;
            ISCONT_AllNetworks(idxObs,7) = 1;
            ISCONT_AllNetworks(idxObs,8) = 1;
            ISCONT_AllNetworks(idxObs,11) = 1;
            ISCONT_AllNetworks(idxObs,12) = 1;
            ISCONT_AllNetworks(idxObs,15) = 1;
            ISCONT_AllNetworks(idxObs,16) = 1;            
        end;    
        for idxStruct = 1:16        
            DeltaISCONT_AllNetworks(idxObs,idxStruct) = (-1)^(ISCONT_AllNetworks(idxObs,idxStruct)+1);
        end;
    end;

    Struct_Indicator_Sum = zeros(1,10);
    for iCase=1:10
        Struct_Indicator_Sum(1,iCase) = sum(Struct_Indicator(:,iCase));
    end;

    %Draw Errors that are used in the calculation of probabilities for all
    %Betas
if (IgnoreStructErrors)
    RealNu2s=Err2;
    Nu2Mean = Mean(RealNu2s)';
    VarNu2=zeros(4,4);
    for n=1:nMkts
        SSample = RealNu2s(n,:)';
        VarNu2 = VarNu2 + ((SSample-Nu2Mean)*(SSample-Nu2Mean)');
    end;
    VarNu2 = VarNu2 / (nMkts);
    CholVarNu2 = chol(VarNu2);
else
    RealNu2s = zeros(nMkts,4);
    Nu2Mean = Mean(RealNu2s)';        
    VarNu2 = zeros(4,4);
    CholVarNu2 = zeros(4,4);
end;    

    EjkMean=mean(Ejk);
    EjkVar=var(Ejk);

    nBetas = (Max-Min)*Grid+1;  

    %======================================================================
    %BEGIN PROBABILITY TESTING
    %======================================================================
    BetaPredictionAll = zeros(nBetaProbIterations,2);
    BetaPredictionAllMeas = zeros(nBetaProbIterations,2);
  if (~doProbMeasErrors)
    for idxIter=1:nBetaProbIterations
        for n=1:nSetsOfErrors
            if ~doProbBootstrap
                Nu2Draw = normrnd(0,1, 4, 1);
                Nu2s(:,n) = Nu2Mean + CholVarNu2' * Nu2Draw;
            else
                RandMktIdx = ceil(rand*nMkts);
                Nu2s(:,n) = RealNu2s(RandMktIdx,:)';
            end;            
        end;

        MomentsPerNetworkStructure = size(MktINST,2);
        ZMoments_AllNetworks = zeros(MomentsPerNetworkStructure,10,nBetas);
        SMoments_AllNetworks = zeros(MomentsPerNetworkStructure,nMkts,10,nBetas);
        [BetaPredictionAll(idxIter,:), ZMoments_AllNetworks, SMomentsAll_AllNetworks, fvalProb]=f_Prob_CalcBeta092605(Min,Max,Grid, A_HMO_All,C_HMO_All,A_Hosp_All,C_Hosp_All,MktINST,ID_Hosp,doVar,ISCONT_AllNetworks,Struct_Indicator,Nu2s,RealNu2s);
        TJTestStatisticProb = fvalProb * nMkts^.5;
        display(BetaPredictionAll);
    end;
    BetaPredMins = sortrows(BetaPredictionAll(:,1));
    BetaPredMaxs = sortrows(BetaPredictionAll(:,2));         

    if doConfInterval
        BetaPredictionInterval(1,1)=BetaPredMins( .025*nBetaProbIterations+1,1 );
        BetaPredictionInterval(2,1)=BetaPredMaxs( .975*nBetaProbIterations,1 );
        display(BetaPredictionInterval);
    else
        display(BetaPredictionAll);
    end;
  end;

    if doProbMeasErrors    
        MeasErrorBetas = zeros(nMeasErrorIter,2);
        MeasErrorSinglePoints=0;

        for idxMeasError=1:nMeasErrorIterProb
            for n=1:nSetsOfErrors
                if ~doProbBootstrap
                    Nu2Draw = normrnd(0,1, 4, 1);
                    Nu2s(:,n) = Nu2Mean + CholVarNu2' * Nu2Draw;
                else
                    RandMktIdx = ceil(rand*nMkts);
                    Nu2s(:,n) = RealNu2s(RandMktIdx,:)';
                end;            
            end;            
            HospMeasurementErrors = normrnd(0,SD_HospMeasurementErrors, nObs/2, 1);
            NBarMeasurementErrors = normrnd(0,SD_NBarMeasurementErrors, nObs/4, 1);
            Etaj = kron(HospMeasurementErrors,ones(2,16));
            Omega = kron(NBarMeasurementErrors,ones(4,16));    

            Psijk_A_All = gamma * (SigmaM_k_All.*SigmaH_jk_All);
            Psijk_B_All = gamma * (SigmaM2_k_All.*Sigma2H_jk_All);
            Psijk_All = Psijk_A_All-Psijk_B_All;
            Psij2k_A_All = gamma * (SigmaM_k_All.*SigmaH_j2k_All);
            Psij2k_B_All = gamma * (SigmaM2_k_All.*Sigma2H_j2k_All);
            Psij2k_All = Psij2k_A_All-Psij2k_B_All;        
            Psijk2_A_All = gamma * (SigmaM_k2_All.*SigmaH_jk2_All);
            Psijk2_B_All = gamma * (SigmaM2_k2_All.*Sigma2H_jk2_All);
            Psijk2_All = Psijk2_A_All-Psijk2_B_All;        

            PsiHosp_A_All = Psijk_A_All + Psijk2_A_All;
            PsiHosp_B_All = Psijk_B_All + Psijk2_B_All;

            DeltaNjkErr_All = DeltaNjk_All + Omega.*Psijk_All;
            AbsDeltaNjkErrInv_All = abs(ones(nObs,16) ./ DeltaNjkErr_All);
            DeltaNjk2Err_All = DeltaNjk2_All + Omega.*Psijk2_All;
            DeltaNj2kErr_All = DeltaNj2k_All + Omega.*Psij2k_All;        

            A_Hosp_All_Meas = DeltaNjkErr_All + DeltaNjk2Err_All;
            A_HMO_All_Meas = (-1) * (DeltaNjkErr_All + DeltaNj2kErr_All);
            C_Hosp_All_Meas = (Njk_All+Njk2_All+Omega.*PsiHosp_A_All).*(HospCostPerPatient_All+Etaj) - (N2jk_All+N2jk2_All+Omega.*PsiHosp_B_All).*(HospCostPerPatient2_All+Etaj) - DeltaNjk_All.*kron(Tjk_Residual,ones(1,16)) - DeltaNjk2_All.*kron(Tjk2_Residual,ones(1,16)) - kron(Delta_Ejk2,ones(1,16));
            C_HMO_All_Meas = (-1)*( ( SigmaM_k_All.*(Prem_All-kron(HMOCost,ones(1,16)))-SigmaM2_k_All.*(Prem2_All-kron(HMOCost,ones(1,16))) ).*(  kron(NBar,ones(1,16)) + Omega) - DeltaNjk_All.*kron(Tjk_Residual,ones(1,16)) - DeltaNj2k_All.*kron(Tj2k_Residual,ones(1,16)) - kron(Delta_Ej2k,ones(1,16)) );
            if ~IgnoreStructErrors
                C_Hosp_All_Meas = C_Hosp_All_Meas - kron(Delta_Ejk,ones(1,16));
                C_HMO_All_Meas = C_HMO_All_Meas - (-1)*(kron(Delta_Ejk,ones(1,16)));
            end;
            MomentsPerNetworkStructureMeas= size(MktINSTmeas,2);
            ZMoments_AllNetworksMeas = zeros(MomentsPerNetworkStructureMeas,10,nBetas);
            SMoments_AllNetworksMeas = zeros(MomentsPerNetworkStructureMeas,nMkts,10,nBetas);
            [BetaPredictionAllMeas(idxMeasError,:), ZMoments_AllNetworksMeas, SMomentsAll_AllNetworksMeas, fvalProbMeas]=f_Prob_CalcBeta092605(Min,Max,Grid, A_HMO_All_Meas,C_HMO_All_Meas,A_Hosp_All_Meas,C_Hosp_All_Meas,MktINSTmeas,ID_Hosp,doVar,ISCONT_AllNetworks,Struct_Indicator,Nu2s,RealNu2s);
            TJTestStatisticProbMeas = fvalProbMeas * nMkts^.5;
            display(BetaPredictionAllMeas);
        end; %foridxMeas        
        BetaPredMeasMins = sortrows(BetaPredictionAllMeas(:,1));
        BetaPredMeasMaxs = sortrows(BetaPredictionAllMeas(:,2));         
        if doConfInterval
            BetaPredictionIntervalMeas(1,1)=BetaPredMeasMins( .025*nMeasErrorIterProb+1,1 );
            BetaPredictionIntervalMeas(2,1)=BetaPredMeasMaxs( .975*nMeasErrorIterProb,1 );
            display(BetaPredictionIntervalMeas);           
        else
            display(BetaPredictionAllMeas);
        end;
    end; %ifUseMeasures


    %======================================================================
    %PROB CONFIDENCE INTERVALS
    %======================================================================
    if doProbVar
        if ~doProbMeasErrors
            BetaConfInterval=zeros(2,1);
            BetaConfInterval2=zeros(2,1);            
            UJSimulAll = zeros(MomentsPerNetworkStructure*10*nBetas,nConfInterval_Iterations);
            ZMomentsAll = reshape(ZMoments_AllNetworks(:,:,:),MomentsPerNetworkStructure*10*nBetas,1);
            VarMoments = zeros(MomentsPerNetworkStructure*10*nBetas);
            for j=1:nMkts
                SSample = [reshape(SMoments_AllNetworks(:,j,:,:),MomentsPerNetworkStructure*10*nBetas,1)];
                VarMoments = VarMoments + ((SSample-ZMomentsAll)*(SSample-ZMomentsAll)')/nMkts;
            end;
            VarMoments = (.99999)*VarMoments + (.00001)*eye(size(VarMoments,1),size(VarMoments,1));            
            BetaSimulation = zeros(nConfInterval_Iterations,2);
            temp = chol(VarMoments);
            for n=1:nConfInterval_Iterations
                epsilons = normrnd(0,1, size(VarMoments,1), 1);
%                UJ_simul = temp' * epsilons / ((nMkts)^.5);
                UJ_simul = temp' * epsilons;
                UJSimulAll(:,n) = UJ_simul;
                Moments_Simul = ZMomentsAll + UJ_simul;
                BetaSimulation(n,:) = f_Prob_CalcBetaFromMoments091305(Moments_Simul,Min,Max,Grid,nObs,MomentsPerNetworkStructure);
            end;
            BetaConfMins = sortrows(BetaSimulation(:,1));
            BetaConfMaxs = sortrows(BetaSimulation(:,2));            
            BetaConfInterval(1,1)=BetaConfMins( .025*nConfInterval_Iterations+1,1 );
            BetaConfInterval(2,1)=BetaConfMaxs( .975*nConfInterval_Iterations,1 );
            display(BetaConfInterval);
       end;
            if doProbMeasErrors    
                BetaConfIntervalMeas=zeros(2,1);
                BetaConfInterval2Meas=zeros(2,1);            
                UJSimulAll = zeros(MomentsPerNetworkStructureMeas*10*nBetas,nConfInterval_Iterations);
                ZMomentsAllMeas = reshape(ZMoments_AllNetworksMeas(:,:,:),MomentsPerNetworkStructureMeas*10*nBetas,1);
                VarMoments = zeros(MomentsPerNetworkStructureMeas*10*nBetas);
                for j=1:nMkts
                    SSample = [reshape(SMoments_AllNetworksMeas(:,j,:,:),MomentsPerNetworkStructureMeas*10*nBetas,1)];
                    VarMoments = VarMoments + ((SSample-ZMomentsAllMeas)*(SSample-ZMomentsAllMeas)')/nMkts;
                end;
                VarMoments = (.99999)*VarMoments + (.00001)*eye(size(VarMoments,1),size(VarMoments,1));            
                BetaSimulation = zeros(nConfInterval_Iterations,2);
                temp = chol(VarMoments);
                for n=1:nConfInterval_Iterations
                    epsilons = normrnd(0,1, size(VarMoments,1), 1);
                    UJ_simul = temp' * epsilons;
                    UJSimulAll(:,n) = UJ_simul;
                    Moments_Simul = ZMomentsAllMeas + UJ_simul;
                    BetaSimulationMeas(n,:) = f_Prob_CalcBetaFromMoments091305(Moments_Simul,Min,Max,Grid,nObs,MomentsPerNetworkStructureMeas);
                end;            
                BetaConfMeasMins = sortrows(BetaSimulationMeas(:,1));
                BetaConfMeasMaxs = sortrows(BetaSimulationMeas(:,2));            
                BetaConfIntervalMeas(1,1)=BetaConfMeasMins( .025*nConfInterval_Iterations+1,1 );
                BetaConfIntervalMeas(2,1)=BetaConfMeasMaxs( .975*nConfInterval_Iterations,1 );
                display(BetaConfIntervalMeas);                
            end; %doConfIntervals for doProbMeasErrors

    end; %doProbVar    
disp(BetaPredictionAllMeas);
end; %Do Probabilities


