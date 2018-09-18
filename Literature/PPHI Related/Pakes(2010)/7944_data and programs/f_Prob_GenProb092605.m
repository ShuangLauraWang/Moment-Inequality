function ProbMatrix=f_Prob_GenProb092605(Beta, A_HMO,C_HMO, A_Hosp,C_Hosp, ISCONTCase1, A_HMO2,C_HMO2, A_Hosp2,C_Hosp2, ISCONTCase2, isMultiStruct, Nu2s,RealNu2s);

    nSetsOfErrors=size(Nu2s,2);
    nObs=size(A_HMO,1);
    nMkts=nObs/4;

    ProbMatrix = zeros(nMkts,1);

    NumStruct = isMultiStruct + 1;
    DeltaH = zeros(nObs,NumStruct);
    DeltaM = zeros(nObs,NumStruct);    
    ISCONT_Both = [ISCONTCase1 ISCONTCase2];
    A_HMO_Both = [A_HMO A_HMO2];
    A_Hosp_Both = [A_Hosp A_Hosp2];
    C_HMO_Both = [C_HMO C_HMO2];
    C_Hosp_Both = [C_Hosp C_Hosp2];

    for idxObs = 1:nObs
        DeltaISCONTCase1(idxObs,1) = (-1)^(ISCONTCase1(idxObs,1)+1);
        DeltaISCONTCase2(idxObs,1) = (-1)^(ISCONTCase2(idxObs,1)+1);        
    end;
    DeltaISCONT_Both = [DeltaISCONTCase1 DeltaISCONTCase2];
    
    for idxStruct = 1:NumStruct    
        DeltaH(:,idxStruct) = (A_Hosp_Both(:,idxStruct) * Beta - C_Hosp_Both(:,idxStruct));  %If Contract, this should be positive
        DeltaM(:,idxStruct) = (A_HMO_Both(:,idxStruct) * Beta - C_HMO_Both(:,idxStruct));    %If Contract, this also should be positive    
    end;
    
    for iErr=1:nSetsOfErrors
        for iMkt=1:nMkts
            ineqfail = zeros(NumStruct,1);            
            ProbTemp = 0;
            iHospAcc = ones(NumStruct,1);
            iHMOAcc = ones(NumStruct,1);
            iHMORej = ones(NumStruct,1);
            
            %Nu2s(:,iErr) = RealNu2s(iMkt,:)';          %In case want to test using real Nu2s
            
            for idxStruct = 1:NumStruct
                for j=1:4                   
                    idx = (iMkt-1)*4+j;

                    if ( (ISCONT_Both(idx,idxStruct)) && ~( DeltaH(idx,idxStruct) + Nu2s(j,iErr) >= 0 ) ) %I(S) - Test if Seller (hosp) accepts, we should have proper inequality
                        iHospAcc(idxStruct,1) = 0;
                        ineqfail(idxStruct,1)=ineqfail(idxStruct,1)+1;
                    end;
                    
                    if ( (ISCONT_Both(idx,idxStruct)) && ~( DeltaM(idx,idxStruct) - Nu2s(j,iErr) >= 0 ) ) %I(B) - Buyer (hmo) accepts
                        iHMOAcc(idxStruct,1) = 0;
                        ineqfail(idxStruct,1)=ineqfail(idxStruct,1)+1;
                    end;

                    if ( ~(ISCONT_Both(idx,idxStruct)) && ~( DeltaM(idx,idxStruct) + Nu2s(j,iErr) >= 0 ) ) %I(B) - Buyer (hmo) rejects
                        iHMORej(idxStruct,1) = 0;
                        ineqfail(idxStruct,1)=ineqfail(idxStruct,1)+1;
                    end;
                end; %for j=1:4
                if ( iHospAcc(idxStruct,1)&&iHMOAcc(idxStruct,1)&&iHMORej(idxStruct,1) )
                    ProbTemp = 1; 
                end;                       
                %display( [Beta iErr iMkt idxStruct ineqfail(idxStruct,1) ProbTemp]);                         
            end; %for NumStruct
            
            ProbMatrix(iMkt,1) = ProbMatrix(iMkt,1) + ProbTemp;
        end; %for nMkts
    end; %for Errs
    nSetsOfErrors = size(Nu2s,2);
    ProbMatrix = ProbMatrix / nSetsOfErrors;