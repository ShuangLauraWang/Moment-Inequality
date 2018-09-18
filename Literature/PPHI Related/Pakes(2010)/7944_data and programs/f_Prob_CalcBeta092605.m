function [BetaPrediction, ZMoments_AllNetworks, SMomentsAll_AllNetworks, fval] = f_Prob_CalcBeta092605(Min,Max,Grid, A_HMO_All,C_HMO_All,A_Hosp_All,C_Hosp_All,INST,ID_Hosp,doVar,ISCONT_AllNetworks,Struct_Indicator,Nu2s,RealNu2s);

nObs = size(A_HMO_All,1);
nMkts=nObs/4;

nBetas = (Max-Min)*Grid + 1;
nInst = size(INST,2);
BetaPrediction = 0;
BetaCurrentMoment = 10^20;

MomentsPerNetworkStructure = nInst;

ZMoments_AllNetworks = zeros(MomentsPerNetworkStructure,10,nBetas);
SMoments_AllNetworks = zeros(MomentsPerNetworkStructure,nMkts,10,nBetas);

BetaMoments=100*ones(nBetas,10);
fval = zeros(nBetas,2);

for idxBeta=1:nBetas
    BetaTemp = (idxBeta-1)/Grid + Min;
    for idxCase=1:10
        [BetaMoments(idxBeta,idxCase), ZMoments_AllNetworks(:,idxCase,idxBeta), SMomentsAll_AllNetworks(:,:,idxCase,idxBeta)]= f_Prob_CalcMoments092605(BetaTemp,idxCase,A_HMO_All,C_HMO_All,A_Hosp_All,C_Hosp_All,INST,doVar,ISCONT_AllNetworks,Struct_Indicator,Nu2s,RealNu2s);
    end; %idxCase=1:10
    %Test what the Value of the sum of moments is across network structures
    TempBetaMomentSum = sum(BetaMoments(idxBeta,:));
    
    display([BetaTemp TempBetaMomentSum]);
    
    if (TempBetaMomentSum < BetaCurrentMoment)
        BetaCurrentMoment = TempBetaMomentSum;
        BetaCurrentPrediction = BetaTemp;          % Set Actual value of Beta to be that of the new idxBeta Min
    else                                            % No longer decreasing in moments; end loop
        %break;
    end;    
end; %idxBeta
fval = BetaCurrentMoment;

if (BetaCurrentMoment == 0)                   %if fval is zero, we don't have unique sol'n. so report all betas that set betamomenttemp=0
    BetaPredictions = [];
    for idxBeta=1:nBetas
        if ( sum(BetaMoments(idxBeta,:)) == 0)
            BetaTemp = (idxBeta-1)/Grid + Min;
            BetaPredictions = [BetaPredictions BetaTemp];
        end;
    end;
    BetaPrediction=[BetaPredictions(1,1) BetaPredictions(1,size(BetaPredictions,2))];
else
    BetaPrediction=[BetaCurrentPrediction BetaCurrentPrediction];
end;

