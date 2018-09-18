function [CalcMoments, ZMomentsStacked, SMomentsStacked] = f_Prob_CalcMoments092605(Beta,idxCase, A_HMO_All,C_HMO_All,A_Hosp_All,C_Hosp_All,INST,doVar,ISCONT_AllNetworks,Struct_Indicator, Nu2s,RealNu2s);

%[BetaMoments(idxBeta,idxCase), ZMoments_AllNetworks(:,idxCase,idxBeta), SMomentsAll_AllNetworks(:,:,idxCase,idxBeta)]
% Return    CalcMoments = Scalar, sum of all negative moments for patricular beta and case
%           ZMomentsStacked = nInst x 1
%           SMomentsStacked = nInst x nMkts

        nObs=size(A_HMO_All,1);    
        nMkts=nObs/4;
        nInst = size(INST,2);
 
        ProbMatrix = zeros(nMkts,1);

        switch idxCase
            case 1
                idxCaseNetwork=1;
                idxCaseNetwork2=1;
            case 2
                idxCaseNetwork=2;
                idxCaseNetwork2=5;
            case 3
                idxCaseNetwork=6;
                idxCaseNetwork2=6;
            case 4
                idxCaseNetwork=3;
                idxCaseNetwork2=9;
            case 5
                idxCaseNetwork=11;
                idxCaseNetwork2=11;
            case 6
                idxCaseNetwork=4;
                idxCaseNetwork2=13;
            case 7
                idxCaseNetwork=7;
                idxCaseNetwork2=10;
            case 8
                idxCaseNetwork=12;
                idxCaseNetwork2=15;
            case 9
                idxCaseNetwork=8;
                idxCaseNetwork2=14;                
            case 10
                idxCaseNetwork=16;
                idxCaseNetwork2=16;                                
        end;
        
        %If current case requires comparison of multiple network
        %structures, e.g. L1S0 has networkstructures 2,5, isMultiStruct=1
        switch idxCase
            case {2,4,6,7,8,9}
                isMultiStruct=1;
            otherwise
                isMultiStruct=0;
        end;
        
        ISCONT_Struct = Struct_Indicator(:,idxCase);
        
        A_Hosp1 = A_Hosp_All(:,idxCaseNetwork);
        C_Hosp1 = C_Hosp_All(:,idxCaseNetwork);
        A_HMO1 = A_HMO_All(:,idxCaseNetwork);
        C_HMO1 = C_HMO_All(:,idxCaseNetwork);
        ISCONTCase1 = ISCONT_AllNetworks(:,idxCaseNetwork);
        A_Hosp2 = A_Hosp_All(:,idxCaseNetwork2);
        C_Hosp2 = C_Hosp_All(:,idxCaseNetwork2);
        A_HMO2 = A_HMO_All(:,idxCaseNetwork2);
        C_HMO2 = C_HMO_All(:,idxCaseNetwork2);
        ISCONTCase2 = ISCONT_AllNetworks(:,idxCaseNetwork2);
        ProbMatrix = f_Prob_GenProb092605(Beta, A_HMO1,C_HMO1, A_Hosp1,C_Hosp1, ISCONTCase1, A_HMO2,C_HMO2, A_Hosp2,C_Hosp2, ISCONTCase2, isMultiStruct, Nu2s,RealNu2s);

        PMoment = (ProbMatrix-ISCONT_Struct);
        ZMoment = zeros(nInst,1);        

        for idxMkt=1:nMkts
                hTemp = INST(idxMkt,:)';
                ZMoment = ZMoment + PMoment(idxMkt,1) * hTemp;                
        end;

        ZMoment = ZMoment/nMkts;        
        SumOfNegativeMoments = 0;
        for i=1:nInst
            if (ZMoment(i,1)<0)
                    SumOfNegativeMoments = SumOfNegativeMoments - ZMoment(i,1);
            end;
        end;

        CalcMoments = SumOfNegativeMoments;
        %display([Beta idxCase CalcMoments]);

        ZMomentsStacked = ZMoment;      %Moment Means
        SMomentsStacked = zeros(nInst,nMkts);                                     %Sample Moments
       
        if doVar
            for iMkt=1:nMkts
                SMoment = zeros(nInst,1);        
                hTemp = INST(idxMkt,:)';
                SMoment = SMoment + PMoment(idxMkt,1) * hTemp;
                SMomentsStacked(:,iMkt) = SMoment;
            end;
        end;