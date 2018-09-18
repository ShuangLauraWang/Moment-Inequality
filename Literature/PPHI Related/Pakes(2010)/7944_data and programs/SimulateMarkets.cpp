#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

// This program simulates 2x2 HMO-Hospital markets and generates
// equilibrium network structures and transfers for randomly drawn
// parameters according to the model described in SimulateMarkets.pdf.
//
// Contact rslee@stern.nyu.edu with any questions.


//=================================================================
//CONSTANTS========================================================
//=================================================================
// Debugging Flags
    const int nRunID = 0;

    const bool bDEBUG = false;                  // Display probability convergence?
    const bool bAllowParameterChanges = false;   // Allow for user to change parameters
    bool bDebugMode=true;                      // Display ALL Transfer Responses?
    bool bIsTransCycling = false;
    bool bIsHMOCycling=false;
    bool bIsNEPrems = true;

    const int nMaxIter = 50000;
    const int nMaxHMOCyclingBeforeQuit=75;     
    const int nMaxCycles=100;                   // 
      
    const double dPercentUseHosp = .075;
    
    const int nHMax = 5;
    const int nMMax = 5;
    const int nSMax = 32;                   // Size of Hospital Strategy Space (2^nMMax))
    const int nRMax = 5;

    int nHospIdx_Large;
    int nHospIdx_Mid;
    int nHospIdx_Small;
    int nIter =1;
    int nSeed =0;

    int nH = 2;                           // Number of Hospitals
    int nM = 2;                           // Number of HMOs
    int nS;
    int nR = 3;
    int nTransUpperBound = 38;
    double dTransGrid = 2;
    int nTransUpperLimitPadding = 0;
    int nTransStratsPerHMO[nHMax];
    int nHMOCyclingBeforeQuit;

// Parameter Distributions
    double nN = 500;                         // Number of Consumers (150 per hosp?)
    double nN_Mean = 600;
    double nN_SD = 300;
    double nN_Min = 100;

    double dCharHosp_Mean = 25;
    double dCharHosp_SD = .5;                //0.5
    double dCostHospConst_Mean = 12;
    double dCostHospConst_SD = 9;
    double dCostHospConst_Min = 1;
    double dCostHMO_Mean = .75;
    double dCostHMO_SD = .25;

    double dCapHosp_Mean = 25;              // 17
    double dCapHosp_SD = 10;                 // 4
    double dCapHosp_Min = 1;
    double dPremium_Mean = 2.25;
    double dPremium_SD = .25;               //.55
    double dCharHMO_Mean = 2.5;
    double dCharHMO_SD = .25;

    double dR_CharHosp_Min = .9;
    double dR_CharHosp_Max = 1.1;
    double dR_CharHMO_Min  = .9;
    double dR_CharHMO_Max  = 1.1;

    double dAlphaPrice = 1.5;
    double dCostHospExpDefault = 1.5;
    double dCorrelationHosp_CharAndCost = .5;
    double dCorrelationHMO_CharAndPrem = .5;

// Global Matricies 
    double dDemCharHMOStd[nMMax];
    int nISCAP[nHMax];

    double dCapHosp[nHMax];
    double dCostHospConst[nHMax];
    double dCostHospExp[nHMax];
    double dCostHMO[nMMax];

    double dTransMax[nMMax];                    // Function of Premium of HMO
    double nTransBounds[nMMax][nHMax][2];       // Lower and Upper Ranges for each HMO-Hospital Transfer pair
    
    double dCharHosp[nHMax];           // Hospital Characteristics
    double dCharHMO[nMMax];            // HMO Characteristics
    double dP[nMMax];                  // Fixed Premiums

    double dR_Shares[nRMax];
    double dR_CharHosp[nRMax][nHMax];
    double dR_CharHMO[nRMax][nMMax];

    double dExpCharHosp[nHMax];
    double dR_ExpCharHosp[nRMax][nHMax];    
    
//=================================================================
//PROTOTYPES=======================================================
//=================================================================
// Main Program Functions
    void vGenerateParameters();
    void vModifyParameters();
        void vPrintParameters();
        void vPrintParameters_File();
        void vChangeParameters();
        double dReadParameter();
    void vCalculateTransfers(double dTrans[nMMax][nHMax]);
    void vPrintDescriptiveStatistics(double dTrans[nMMax][nHMax]);
    void vPrintDescriptiveStatistics_File(double dTrans[nMMax][nHMax]);
    void vOutputData(double dTrans[nMMax][nHMax]);
    void vOutputData2(int nOutputType, int nHospIdx, double dTrans[nMMax][nHMax], double dPrem[nMMax], int nHMOStruct[nMMax], int nHospStruct[nHMax]);
    void vOutputData3(int nHospIdx, double dTrans[nMMax][nHMax], int nHMOStruct[nMMax],int nHospStruct[nHMax]);
    void vCheckIsNEPrems(double dTrans[nMMax][nHMax]);
    void vCopyStruct(int nStructCopy[nHMax], int nHospStruct[nHMax]);

    void vPrintHospitalNormalForm(double dTrans[nMMax][nHMax]);

// Stage 4
    double dSigmaM(int nMIndex, int nHospStruct[nHMax], double dPrem[nMMax]);
    double dSigmaHM(int nHidx, int nMidx, int nHospStruct[nHMax], double dPrem[nMMax]);
    double dSigmaH_R (int nHidx, int nRidx, int nMidx, int nHospStruct[nHMax],double dPrem[nMMax]);
    double dSigmaM_R (int nMidx, int nRidx, int nHospStruct[nHMax],double dPrem[nMMax]);
    double dSigmaR_M (int nMidx, int nRidx, int nHospStruct[nHMax],double dPrem[nMMax]);

//Stage 3
    void vSetPremiums(double dPrem[nMMax], int nHospStruct[nHMax], double dTrans[nMMax][nHMax]);
    double dPremiumFOC(int nMidx, double dPrem[nMMax], int nHospStruct[nHMax], double dTrans[nMMax][nHMax]);
    double dPHosp (int nHidx, double dTrans[nMMax][nHMax]);
    double dPHosp_Struct (int nHidx, int nHospStruct[nMMax], double dTrans[nMMax][nHMax]);
    double dCostsHosp(int nHidx, int nHospStruct[nMMax], double dPrem[nMMax]);
    double dPHMO(int nMidx, int nHMOStruct[nMMax], double dTrans[nMMax][nHMax]);
    double dPHMO_HospStruct(int nMidx, int nHospStruct[nHMax], double dTrans[nMMax][nHMax]);
    void vHMO_GenNE(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]);
    void vHMO_GenNE_fullsearch(int nHMOStruct[nHMax],double dTrans[nMMax][nHMax]);
    bool bHMO_IsNE(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]);
    void vHMO_GenBR(int nHMOStruct[nHMax],double dTrans[nMMax][nHMax]);
    void vConvertStruct_HMO2Hosp(int nHMOStruct[nHMax], int nHospStruct[nMMax]);
    void vChangeStruct_HM(int nHidx, int nMidx, int nHospStruct[nHMax], int nHospStructNew[nHMax]);
    void vChangeStruct_HM2(int nHidx, int nMidx, int nHospStruct[nHMax], int nHospStructNew[nHMax]);
    
    void vCopyTransfers(double dTransCopy[nMMax][nHMax], double dTrans[nMMax][nHMax]);
    void vPrintTransfers(double dTrans[nMMax][nHMax]);
    void vPrintTransfers_File(double dTrans[nMMax][nHMax]);
    void vPrintHMOStructure_FromTrans(double dTrans[nMMax][nHMax]);
    void vPrintHMOStructure(int nHMOStruct[nMMax]);
    void vPrintHospStructure_FromTrans(double dTrans[nMMax][nHMax]);
    void vPrintHospStructure(int nHospStruct[nHMax]);

    int nIndustryStructure_FromHMOStruct(int nHMOStruct[nMMax]);
    void vPrintIndustryStructure_FromHMOStruct(int nHMOStruct[nMMax]);
    void vPrintIndustryStructure_FromHMOStruct_File(int nHMOStruct[nMMax]);

//Probabilities
    double normal_01_cdf ( double x );
    double normal_01_cdf_inv ( double cdf );
    double d_huge ( void );
    double dpoly_value ( int n, double a[], double x );
    double normal_01_sample ( int *seed );
    int get_seed ( void );
    double d_uniform_01 ( int *seed );
    
//=================================================================
//MAIN PROGRAM=====================================================
//=================================================================
int main()
{
    while (nIter < nMaxIter) {
        cout << "ITERATION: " << nIter << "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
        ofstream oFile("output.txt",ios::app);
        oFile << "ITERATION: " << nIter << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
        oFile.close();

        nSeed = get_seed();
        
        vGenerateParameters();
    
        vPrintParameters();
        vPrintParameters_File();    

        if (bAllowParameterChanges) vModifyParameters();
    
        for (int i=0; i<nHMax; i++) {
            dExpCharHosp[i] = exp(dCharHosp[i]);           //Calculate Exp(Characteristics) for computational speed
            for (int r=0; r<nR; r++) dR_ExpCharHosp[r][i] = exp(dR_CharHosp[r][i] * dCharHosp[i]);
        }
        
        //======StartingTransfers========================
        cout << "Transfer Ranges for Each Hospital:\n";
        double dTrans[nMMax][nHMax];
        for (int j=0;j<nH;j++) {
           for (int k=0;k<nM;k++) {
                dTrans[k][j]=0;
                dTrans[k][j]= (int(dCostHospConst[j]/dTransGrid) + 1)*(dTransGrid);
            }
            nTransStratsPerHMO[j] = (int) (nTransUpperBound / dTransGrid);
            cout << "[" << nTransBounds[0][j][0] << "," << (nTransBounds[0][j][0]+nTransStratsPerHMO[j]-1)*dTransGrid << "]\t";
            
        }
        cout << endl;
            
        //=======Main Program=============================== 
        int start=time(NULL);
        vCalculateTransfers(dTrans);
        if (!bIsHMOCycling&&!bIsTransCycling) {
            vPrintDescriptiveStatistics(dTrans);
            vPrintDescriptiveStatistics_File(dTrans);
            vCheckIsNEPrems(dTrans);
            vOutputData(dTrans);
        }
        oFile.open("output.txt",ios::app);
        oFile << "Runtime: " << time(NULL) - start << " seconds.\n\n";
        oFile.close();
        cout << "Runtime: " << time(NULL) - start << " seconds.\n\n";
        nIter++;
    }
    return 0;
}


//=================================================================
//FUNCTIONS=====================================================
//=================================================================

void vGenerateParameters() {
    bIsTransCycling=false;
    bIsHMOCycling=false;
    nHMOCyclingBeforeQuit = 0;
    nS = (int)pow(2.,nH);
    double dRho = dCorrelationHosp_CharAndCost;

    nN = nN_Mean + normal_01_sample(&nSeed) * nN_SD;
    if (nN < nN_Min) nN = nN_Min;

    //Hospital Parameters;
    for (int i=0; i<nHMax; i++) {
        dCharHosp[i] = normal_01_sample(&nSeed);
        dCostHospConst[i] = dCostHospConst_SD*(dRho*dCharHosp[i] + pow((1-pow(dRho,2.)),.5) * normal_01_sample(&nSeed)) + dCostHospConst_Mean;
        if (dCostHospConst[i] < dCostHospConst_Min) dCostHospConst[i] = dCostHospConst_Min;
        dCharHosp[i] = dCharHosp[i] * dCharHosp_SD + dCharHosp_Mean;
        dCostHospExp[i] = dCostHospExpDefault;
        dCapHosp[i] = normal_01_sample(&nSeed) * dCapHosp_SD + dCapHosp_Mean;
        if (dCapHosp[i] < dCapHosp_Min) dCapHosp[i]=dCapHosp_Min;
    }
    double dCapTemp;
    if (dCapHosp[0] < dCapHosp[1]) {
        dCapTemp = dCapHosp[0];
        dCapHosp[0] = dCapHosp[1];
        dCapHosp[1] = dCapTemp;
    }
    
    //HMO Parameters
    dRho = dCorrelationHMO_CharAndPrem;
    for (int i=0; i<nMMax; i++) {
        dP[i] = normal_01_sample(&nSeed);
        dCharHMO[i] = dCharHMO_SD*(dRho*dP[i] + pow((1-pow(dRho,2.)),.5) * normal_01_sample(&nSeed)) + dCharHMO_Mean;
        dCostHMO[i] = dP[i] * dCostHMO_SD + dCostHMO_Mean;
        dP[i] = dP[i] * dPremium_SD + dPremium_Mean;
        dTransMax[i] = dP[i] / dPercentUseHosp;
    }
    //Demographic Parameters
    double dR_SharesSum = 0;
    for (int r=0; r<nR; r++) {
        for (int j=0; j<nH; j++) {
            dR_CharHosp[r][j] = d_uniform_01(&nSeed)*(dR_CharHosp_Max-dR_CharHosp_Min) + dR_CharHosp_Min;
        }
        for (int k=0; k<nM; k++) {
            dR_CharHMO[r][k] = d_uniform_01(&nSeed)*(dR_CharHMO_Max-dR_CharHMO_Min) + dR_CharHMO_Min;;
        }
        dR_Shares[r] = d_uniform_01(&nSeed);
        dR_SharesSum = dR_SharesSum + dR_Shares[r];
    }
    for (int r=0; r<nR; r++) dR_Shares[r] = dR_Shares[r] / dR_SharesSum;   
    for (int k=0; k<nM; k++) {
        dDemCharHMOStd[k] = 0;
        for (int r=0; r<nR; r++) {
            dDemCharHMOStd[k] = dDemCharHMOStd[k] + dR_Shares[r]*dR_CharHMO[r][k]*dR_CharHMO[r][k];
        }
        double dTemp = 0;
        for (int r=0; r<nR; r++) {
            dTemp = dTemp + dR_Shares[r]*dR_CharHMO[r][k];        
        }
        dDemCharHMOStd[k] = dDemCharHMOStd[k] - dTemp*dTemp;
        dDemCharHMOStd[k] = pow(dDemCharHMOStd[k],.5);
    }     
    
    //Set Bounds on Transfer Space
    int nTemp  = 0;
    int nTemp2 = 0;
    for (int j=0;j<nH;j++) {
        nTemp = 0;
        for (int k=0; k<nM; k++) {
//            nTransBounds[k][j][0] = (int)(dCostHospConst[j]/dTransGrid) - nTransLowerLimitPadding;
//            if (nTransBound[k][j][0] < 0) nTransBound[k][j][0] =0;
            nTransBounds[k][j][0] = 0;
            nTransBounds[k][j][1] = (int)(dTransMax[k]/dTransGrid) + nTransUpperLimitPadding;
            nTemp2 = (int) (nTransBounds[k][j][1] - nTransBounds[k][j][0] + 1);
            if (nTemp2 > nTemp) nTemp = nTemp2;
        }
        nTransStratsPerHMO[j] = nTemp;
    }
    cout << endl;
    
    //HospitalSizes
    nHospIdx_Large=0;
    for (int j=0; j<nH; j++) {
        if (dCapHosp[j] > dCapHosp[nHospIdx_Large]) nHospIdx_Large=j;        
    }
    nHospIdx_Mid=(nHospIdx_Large+1)%nH;
    for (int j=1; j<nH; j++) {
        if (dCapHosp[(nHospIdx_Large+j)%nH] > dCapHosp[(nHospIdx_Mid)%nH]) nHospIdx_Mid=(nHospIdx_Large+j)%nH;
    }    
    nHospIdx_Small=0;
    for (int j=0; j<nH; j++) {
        if (nHospIdx_Small==nHospIdx_Large) nHospIdx_Small++;
        if (nHospIdx_Small==nHospIdx_Mid) nHospIdx_Small++;
    }
    cout << nHospIdx_Large << " " << nHospIdx_Mid << " " << nHospIdx_Small << endl;
}

void vModifyParameters() {
    if (bAllowParameterChanges) {
        bool bChangeParameters = true;
        int nInput=0;
        while (bChangeParameters) {
            bChangeParameters = false;
            cout << "Change Parameters? -- Enter 1 for yes, 0 for no: ";
            cin >> nInput;
            cin.ignore();
            if (nInput) {
                vChangeParameters();
                vPrintParameters();
                bChangeParameters = true;
            }
        }
        cout << "Debugging Mode (Full Transfer Display)? -- Enter 1 for yes, 0 for no: ";
        cin >> nInput;
        cin.ignore
        ();
        if (nInput) {
            bDebugMode = true;
        }
    }
}

void vPrintParameters() {
    cout.precision(4);
    cout << "======================" << endl;
    cout << "Parameters of Model" << endl;
    cout << "======================" << endl;
    cout << "Pop =" << nN << "\t #HMO =" << nM <<"\t #HOS =" << nH << endl;
    for (int i=0;i<nH;i++) {
        cout << "HOS "<< i+1 << ": \t Char =" << dCharHosp[i] << "\t Costs =" << dCostHospConst[i] << " , " << dCostHospExp[i] << "\t Capacity =" << dCapHosp[i] << endl;
    }
    for (int i=0;i<nM;i++) {
        cout << "HMO "<< i+1 << ": \t Char =" << dCharHMO[i] << "\t Costs =" << dCostHMO[i] << "\t Net =" << dCharHMO[i]-dAlphaPrice*dP[i] << endl;
    }
    cout << "Transfer Increments: "<< dTransGrid << endl;
    cout << endl;
    cout << "Demographic Shares: \t";
    for (int r=0; r<nR; r++) {
        cout << dR_Shares[r] << "\t";
    }
    cout << endl;
    cout << "Demographic Prefs: \t H1\t H2\t H3\t M1\t M2\t M3" << endl;
    for (int r=0; r<nR; r++) {
        cout << "                  \t ";
        for (int i=0; i<nH; i++) cout << dR_CharHosp[r][i] << "\t";
        for (int i=0; i<nM; i++) cout << dR_CharHMO[r][i] << "\t";
        cout << endl;
    }
    cout << endl;
    cout << "StdDev in HMO Char: \t";
    for (int i=0; i<nM; i++) cout << dDemCharHMOStd[i] << "\t";
    cout << endl;
}

void vPrintParameters_File() {
    ofstream oFile("output.txt",ios::app);
    oFile.precision(4);
    oFile << "======================" << endl;
    oFile << "Parameters of Model" << endl;
    oFile << "======================" << endl;
    oFile << "Pop =" << nN << "\t #HMO =" << nM <<"\t #HOS =" << nH << endl;
    for (int i=0;i<nH;i++) {
        oFile << "HOS "<< i+1 << ": \t Char =" << dCharHosp[i] << "\t Costs =" << dCostHospConst[i] << " , " << dCostHospExp[i] << "\t Capacity =" << dCapHosp[i] << endl;
    }
    for (int i=0;i<nM;i++) {
        oFile << "HMO "<< i+1 << ": \t Char =" << dCharHMO[i] << "\t Costs =" << dCostHMO[i] << "\t Net =" << dCharHMO[i]-dP[i] << endl;
    }
    oFile << "Transfer Increments: "<< dTransGrid << endl;
    oFile << endl;
    oFile << "Demographic Shares: \t";
    for (int r=0; r<nR; r++) {
        oFile << dR_Shares[r] << "\t";
    }
    oFile << endl;
    oFile << "Demographic Prefs: \t H1\t H2\t H3\t M1\t M2\t M3" << endl;
    for (int r=0; r<nR; r++) {
        oFile << "                  \t ";
        for (int i=0; i<nH; i++) oFile << dR_CharHosp[r][i] << "\t";
        for (int i=0; i<nM; i++) oFile << dR_CharHMO[r][i] << "\t";
        oFile << endl;
    }
    oFile << endl;
    oFile << "StdDev in HMO Char: \t";
    for (int i=0; i<nM; i++) oFile << dDemCharHMOStd[i] << "\t";
    oFile << endl;
    oFile.close();
}

void vChangeParameters() {
    double dTemp;

    // Set N, M, H
        cout << "Population (N) (default=" << nN <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) nN = (int)dTemp;

        cout << "# HMOs (M) (default=" << nM <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) nM = (int)dTemp;

        cout << "# Hospitals (H) (default=" << nH <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) nH = (int)dTemp;

        nS = (int)pow(2.,nM);

    // Set Others
    for (int i=0;i<nH;i++) {
        cout << "HOS " << i+1 << " Characteristic (default=" << dCharHosp[i] <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) dCharHosp[i] = dTemp;
    }
    for (int i=0;i<nH;i++) {
        cout << "HOS " << i+1 << " Cost Per Patient (default=" << dCostHospConst[i] <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) dCostHospConst[i] = dTemp;
    }
    for (int i=0;i<nH;i++) {
        cout << "HOS " << i+1 << " Capacity (default=" << dCapHosp[i] <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) dCapHosp[i] = dTemp;
    }
    for (int i=0;i<nM;i++) {
        cout << "HMO " << i+1 << " Characteristic (default=" << dCharHMO[i] <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) dCharHMO[i] = dTemp;
    }
    for (int i=0;i<nM;i++) {
        cout << "HMO " << i+1 << " Fixed Premium (default=" << dP[i] <<"): ";
        dTemp=dReadParameter();
        if (dTemp>0) dP[i] = dTemp;
    }

    cout << "Transfers Unit (default=" << dTransGrid << "): ";
    dTemp=dReadParameter();
    if (dTemp>0) dTransGrid = dTemp;
}

double dReadParameter() {
    double dTemp = -1;
    char szTemp[50];
    while (dTemp<0) {
        cin.getline(szTemp,50);
        if (!(*szTemp)) return -1;
            else dTemp=atof(szTemp);
    }
    return dTemp;
}

void vCalculateTransfers(double dTrans[nMMax][nHMax]) {
    ofstream oFile("output.txt",ios::app);
    oFile << "Now Calculating Transfers" << endl;
    oFile.close();

    cout << "======================" << endl;
    cout << "Now Calculating Transfers" << endl;
    cout << "======================" << endl;
    double dTransTemp[nMMax][nHMax];
    double dTransHistory[nMaxCycles][nMMax][nHMax];
    int nCycles=0;
    bool bIsFP = false;
    vPrintTransfers(dTrans);
    vPrintTransfers_File(dTrans);    
    while (!bIsFP) {
        bIsFP=true;
        
        for (int j=0;j<nH;j++) {                // look at all Hosps
            double dTMax = (nTransStratsPerHMO[j]-1) * dTransGrid; 
            double dCurrentHospProfits = dPHosp(j, dTrans);
            double dNewHospProfits = dCurrentHospProfits;
            vCopyTransfers(dTransTemp,dTrans);
            int nSizeOfTransferSpacePerHMO = nTransStratsPerHMO[j];           
            for (int nTransState=0; nTransState<(int)pow((double)nSizeOfTransferSpacePerHMO,nM); nTransState++) {
                bIsHMOCycling = false;
                for (int k=0;k<nM;k++) {
                    dTransTemp[k][j] = (((nTransState % (int)pow((double)nSizeOfTransferSpacePerHMO,nM-k)) / (int)pow((double)nSizeOfTransferSpacePerHMO,nM-1-k)) + (double)nTransBounds[k][j][0] ) * dTransGrid;
                }
                dNewHospProfits=dPHosp(j, dTransTemp);
                if (bIsHMOCycling) {
                    dNewHospProfits = 0;
                    nHMOCyclingBeforeQuit++;
                    if (nHMOCyclingBeforeQuit > nMaxHMOCyclingBeforeQuit) goto STOP;
                }
                
                if (dNewHospProfits > dCurrentHospProfits) {
                    bIsFP=false;
                    dCurrentHospProfits=dNewHospProfits;
                    for (int k=0;k<nM;k++) {
                        dTrans[k][j] = dTransTemp[k][j];
                    }
                }
            }
            for (int k=0;k<nM;k++) {
                if (dTrans[k][j]==dTMax) nTransStratsPerHMO[j]=nTransStratsPerHMO[j]+5;
            }  
            if (bDebugMode) {
                cout << "\t";
                vPrintTransfers(dTrans);
                oFile.open("output.txt",ios::app);
                oFile << "\t";
                oFile.close();
                vPrintTransfers_File(dTrans);
                vPrintHospStructure_FromTrans(dTrans);
            }
        }
        vPrintTransfers(dTrans);
        vPrintTransfers_File(dTrans);
        int nTransCycleLength;
        if (nCycles==nMaxCycles) {
            cout << "No Convergence in Transfers!" << endl;
            break;
        } else {
            vCopyTransfers(dTransHistory[nCycles],dTrans);
            for (int i=0;i<nCycles-1;i++) {
                bIsTransCycling=true;
                for (int k=0; k<nM; k++) {
                    for (int j=0; j<nH; j++) {
                        if (dTransHistory[i][k][j]!=dTrans[k][j]) {
                            bIsTransCycling=false;
                        }
                    }
                }
                if (bIsTransCycling) {
                    nTransCycleLength = nCycles - i;
                    break;
                }
            }
        }
        if (bIsTransCycling) {
            cout << "ERROR: Transfers are cycling with length: " << nTransCycleLength << endl;
            oFile.open("output.txt",ios::app);
            oFile << "ERROR: Transfers are cycling with length: " << nTransCycleLength << endl;
            oFile.close();
            break;
        }
        nCycles++;
    }
    STOP: if (bIsHMOCycling) cout << "Procedure Terminated";
}

void vPrintDescriptiveStatistics (double dTrans[nMMax][nHMax]) {
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);    
    double dPrem[nMMax];
    vSetPremiums(dPrem,nHospStruct,dTrans);

    cout << "======================" << endl;
    cout << "Descriptive Statistics" << endl;
    cout << "======================" << endl;

    cout << "Final Trans :" << endl;
    vPrintTransfers(dTrans);
 
    cout.precision(2);
    cout << "Premiums :" << dPrem[0] << " " << dPrem[1] << endl;
    
    cout << "Realized Profits:" << endl;
    cout.precision(3);
    cout << "\t HMO : \t";    
    for (int k=0;k<nM;k++) {
        cout << dPHMO(k, nHMOStruct, dTrans) << "\t";
    }
    cout << endl;
    cout << "\t HOS : \t";    
    for (int i=0;i<nH;i++) {
        cout << dPHosp(i, dTrans) << "\t";\
    }
    cout << endl;

    vPrintIndustryStructure_FromHMOStruct(nHMOStruct);

    cout << "HMO Shares : \t";    
    cout.precision(3);
    for (int k=0;k<nM;k++) cout << dSigmaM(k, nHospStruct, dPrem) << "\t";
    cout << endl << endl;    

    cout << "Hosp Shares : \n";    
    cout.precision(3);
    for (int j=0;j<nH;j++) {
        cout <<"\t Hosp" << j+1 << ": \t";
        for (int k=0;k<nM;k++) {
            cout << dSigmaHM(j, k, nHospStruct, dPrem) / dPercentUseHosp << "\t";        
        }
        cout << endl;
    }
    cout << endl << endl;    
    
    cout << "Transfer Cycling :" << bIsTransCycling << endl;
    cout << "HMO Cycling :" << bIsHMOCycling << endl;
    cout << endl;
}

void vPrintDescriptiveStatistics_File (double dTrans[nMMax][nHMax]) {
    ofstream oFile("output.txt",ios::app);
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);
    double dPrem[nMMax];
    vSetPremiums(dPrem,nHospStruct,dTrans);

    oFile << "======================" << endl;
    oFile << "Descriptive Statistics" << endl;
    oFile << "======================" << endl;

    oFile << "Final Trans :" << endl;
    oFile.close();
    vPrintTransfers_File(dTrans);
    oFile.open("output.txt",ios::app); 

    oFile << "Premiums : \t";
    for (int k=0;k<nM;k++) {
        oFile << dPrem[k] << "\t";
    }
    oFile << endl;
    
    oFile << "Realized Profits:" << endl;
    oFile << "\t HMO : \t";    
    for (int k=0;k<nM;k++) {
        oFile << dPHMO(k, nHMOStruct, dTrans) << "\t";
    }
    oFile << endl;
    oFile << "\t HOS : \t";    
    for (int i=0;i<nH;i++) {
        oFile << dPHosp(i, dTrans) << "\t";\
    }
    oFile << endl;
    oFile.close();
    vPrintIndustryStructure_FromHMOStruct_File(nHMOStruct);
    oFile.open("output.txt",ios::app);
    oFile << "HMO Shares : \t";    
    for (int k=0;k<nM;k++) oFile << dSigmaM(k, nHospStruct, dPrem) << "\t";
    oFile << endl << endl;
    oFile << "Hosp Shares : \n";    
    for (int j=0;j<nH;j++) {
        oFile <<"\t Hosp" << j+1 << ": \t";
        for (int k=0;k<nM;k++) {
            oFile << dSigmaHM(j, k, nHospStruct, dPrem)/dPercentUseHosp << "\t";
        }
        oFile << endl;
    }
    oFile << endl << endl;    
    
    oFile << "Transfer Cycling :" << bIsTransCycling << endl;
    oFile << "HMO Cycling :" << bIsHMOCycling << endl;
    oFile << endl;
    oFile.close();
}

void vOutputData(double dTrans[nMMax][nHMax]) {
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);
    double dPrem[nMMax];
    vSetPremiums(dPrem,nHospStruct,dTrans);

    //ISCAP Matrix
    for (int j=0; j<nH;j++) {
        //Tests to see if hospital j is at capacity when market is fully contracted
        int nHospStruct2[nHMax];
        double dNMH2 = 0;
        for (int l=0; l<nHMax;l++) nHospStruct2[l]=nS-1;
        for (int k=0;k<nM;k++) {    
            double dHMOShare2 = dSigmaM(k,nHospStruct2,dPrem);
            double dHospShare2 = dSigmaHM(j,k,nHospStruct2,dPrem);
            dNMH2 = dNMH2 + nN * dHMOShare2 * dHospShare2;
        }
        if (dNMH2 > dCapHosp[j]) {
            nISCAP[j] = 1;
        } else {
            nISCAP[j] = 0;
        }        
    }
    
    vOutputData2(0, nHospIdx_Large, dTrans, dPrem, nHMOStruct, nHospStruct);
    vOutputData2(0, nHospIdx_Mid, dTrans, dPrem, nHMOStruct, nHospStruct);
    vOutputData3(nHospIdx_Large, dTrans, nHMOStruct,nHospStruct);
    vOutputData3(nHospIdx_Mid, dTrans, nHMOStruct,nHospStruct);    
}

void vOutputData2(int nOutputType, int nHospIdx, double dTrans[nMMax][nHMax], double dPrem[nMMax], int nHMOStruct[nMMax], int nHospStruct[nHMax]) {
    char strOutputFile[30];
    switch ( nOutputType )
    {
        case 0: 
            strcpy(strOutputFile,"data.txt");
            break;
        case 1: 
            strcpy(strOutputFile,"data_HospA.txt");
            break;
        case 2: 
            strcpy(strOutputFile,"data_HospB.txt");
            break;
        case 3: 
            strcpy(strOutputFile,"data_HospC.txt");
            break;
    }
    ofstream oFile(strOutputFile,ios::app);

    //OUTPUT DATA==========================================================
    int j = nHospIdx;
    
        double dRealizedHospProfits = dPHosp_Struct(j,nHospStruct,dTrans);
        double dRealizedHospCosts = dCostsHosp(j,nHospStruct,dPrem);

        for (int k=0;k<nM;k++) {
            oFile << nIter << "\t" << k << "\t" << j << "\t" << nN << "\t";        // [ITER#] , [HMO_ID], [HOSP_ID], MktSize

            //HMOs JOINED BY HOSPITAL                                                               // ISCONT
            if (nHospStruct[j]&(int)pow(2.,k)) {
                oFile << "1 \t";
            } else {
                oFile << "0 \t";
            }
            
            //AllHospCosts                                                                          // [HospCosts: C_j C_j' C_j'']
            for (int nTemp=0; nTemp<nH; nTemp++) {
                oFile << dCostHospConst[(j+nTemp)%nH] << "\t";
            }            

            //AllISCAP                                                                              // [ISCAP: I_j I_j' I_j'']
            for (int nTemp=0; nTemp<nH; nTemp++) {
                oFile << nISCAP[(j+nTemp)%nH] << "\t";
            }      

            //AllHospCap                                                                            // [HospCap: Gam_j Gam_j' Gam_j'']
            for (int nTemp=0; nTemp<nH; nTemp++) {
                oFile << dCapHosp[(j+nTemp)%nH] << "\t";
            }

            oFile << dCharHosp[j] << "\t";                                                          // [HOSPCHAR]
    

            //TRANSFERS OFFERED BY HOSPITAL                                                         // [TRANSFERS: T_jk T_jk' T_jk'']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                oFile << dTrans[(k+nTemp)%nM][j] << "\t";
            }
            for (int nTemp=1; nTemp<nH; nTemp++) {                                                  // [TRANSFERS: T_j'k T_j''k]
                oFile << dTrans[k][(j+nTemp)%nH] << "\t";
            }
            
            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO k,k',k''                                // [N_jk, N_jk', N_jk'']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                int nTempHMOidx = (k+nTemp)%nM;
                double dHMOShare = dSigmaM(nTempHMOidx,nHospStruct,dPrem);
                double dHospShare = dSigmaHM(j,nTempHMOidx,nHospStruct,dPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            
            for (int nTemp=1; nTemp<nH; nTemp++) {                                                  // [N_j'k, N_j''k]
                int nTempHospidx = (j+nTemp)%nH;
                double dHMOShare = dSigmaM(k,nHospStruct,dPrem);
                double dHospShare = dSigmaHM(nTempHospidx,k,nHospStruct,dPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            
            
            //Realized Hospital Profits and Costs                                                   // [RealizedHospProfits, RealizedHospCosts]
            oFile << dRealizedHospProfits << "\t" << dRealizedHospCosts << "\t";
            
            //Alternative Structure -- Hospital Profits and Costs and Prems                                  // [NewHospProfits,NewHospCosts,NewPrem]
            int nHospStruct_New[nHMax];
            vChangeStruct_HM(j,k,nHospStruct,nHospStruct_New);
            double dNewPrem[nMMax];
            vSetPremiums(dNewPrem,nHospStruct_New,dTrans);
            oFile << dPHosp_Struct(j,nHospStruct_New,dTrans) << "\t";
            oFile << dCostsHosp(j,nHospStruct_New,dNewPrem) << "\t";
            oFile << dNewPrem[k] << "\t";

            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO k,k',k'', Alternative                   // [N'_jk, N'_jk', N'_jk'']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                int nTempHMOidx = (k+nTemp)%nM;
                double dHMOShare = dSigmaM(nTempHMOidx,nHospStruct_New,dNewPrem);
                double dHospShare = dSigmaHM(j,nTempHMOidx,nHospStruct_New,dNewPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            

            //HMO Stuff ========================================================================
            oFile << dCharHMO[k] << "\t";                                                           // [HMOChar]
            oFile << dPrem[k] << "\t";                                                              // [HMOPremiums]
            oFile << dCostHMO[k] << "\t";                                                           // [HMO Cost]
            oFile << dPHMO_HospStruct(k,nHospStruct,dTrans) << "\t";                                // [HMOProfits]

            for (int nTemp=0; nTemp<nM; nTemp++) {                                                  // [SigmaM_k, SigmaM_k', SigmaM_k'']
                oFile << dSigmaM((k+nTemp)%nM, nHospStruct, dPrem) << "\t";
            }
            
            oFile << dPHMO_HospStruct(k,nHospStruct_New,dTrans) << "\t";                            // [NewHMOProfits]
            for (int nTemp=0; nTemp<nM; nTemp++) {                                                  // [SigmaM'_k, SigmaM'_k', SigmaM'_k'']
                oFile << dSigmaM((k+nTemp)%nM, nHospStruct_New, dNewPrem) << "\t";
            }
            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO                                         // [N'_jk, N'_j'k, N'_j''k]
            for (int nTemp=0; nTemp<nH; nTemp++) {
                int nTempHospidx = (j+nTemp)%nH;
                double dHMOShare = dSigmaM(k,nHospStruct_New,dNewPrem);
                double dHospShare = dSigmaHM(nTempHospidx,k,nHospStruct_New,dNewPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            
            //DemStdDevHMOChar                                                                      // [TRANSFERS: DemHMOCharSTD_k DemHMOCharSTD_k' DemHMOCharSTD_k'']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                oFile << dDemCharHMOStd[(k+nTemp)%nM] << "\t";
            }      
            oFile << endl;
        }        
        
    oFile.close();
}

void vOutputData3(int nHospIdx, double dTrans[nMMax][nHMax], int nHMOStruct[nMMax],int nHospStruct[nMMax]) {
    char strOutputFile[30];
    strcpy(strOutputFile,"data_AllNetworkStructures.txt");
    ofstream oFile(strOutputFile,ios::app);
    int j = nHospIdx;
    for (int k=0;k<nM;k++) {    

        //OUTPUT DATA==========================================================
        oFile << nIter << "\t" << k << "\t" << j << "\t";                                       // [ITER#] , [HMO_ID], [HOSP_ID]
        oFile << nIndustryStructure_FromHMOStruct(nHMOStruct) << "\t";                          // [EQUILIBRIUM NETWORK STRUCTURE]

        //TRANSFERS OFFERED BY HOSPITAL                                                         // [TRANSFERS: T_jk T_jk']
        for (int nTemp=0; nTemp<nM; nTemp++) {
            oFile << dTrans[(k+nTemp)%nM][j] << "\t";
        }
        for (int nTemp=1; nTemp<nH; nTemp++) {                                                  // [TRANSFERS: T_j'k]
            oFile << dTrans[k][(j+nTemp)%nH] << "\t";
        }
       
        for (int nC = 0; nC < (int)pow((double)nS,nM); nC++) {                                  // REPEAT FOR ALL NETWORK STRUCTURES:::
            oFile << nC << "\t";                                                                // NETWORKSTRUCTURE_ID (0-15)
            int nHMOStructAlt[nMMax];
            for (int i=0; i<nM; i++) {
                nHMOStructAlt[i] = (nC % (int)pow((double)nS,nM-i)) / (int)pow((double)nS,nM-1-i);
            }
            vConvertStruct_HMO2Hosp(nHMOStructAlt,nHospStruct);
            double dPrem[nMMax];
            vSetPremiums(dPrem,nHospStruct,dTrans);
            int nHospStruct_New[nHMax];
            vChangeStruct_HM(j,k,nHospStruct,nHospStruct_New);
            double dNewPrem[nMMax];
            vSetPremiums(dNewPrem,nHospStruct_New,dTrans);

            double dRealizedHospProfits = dPHosp_Struct(j,nHospStruct,dTrans);
            double dRealizedHospCosts = dCostsHosp(j,nHospStruct,dPrem);
            double dNewHospProfits =   dPHosp_Struct(j,nHospStruct_New,dTrans);
            double dNewHospCosts  = dCostsHosp(j,nHospStruct_New,dNewPrem);

            //HMOs JOINED BY HOSPITAL                                                               // ISCONT[j,k]
            if (nHospStruct[j]&(int)pow(2.,k)) {
                oFile << "1 \t";
            } else {
                oFile << "0 \t";
            }
        
            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO k,k',k''                                // [N_jk, N_jk']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                int nTempHMOidx = (k+nTemp)%nM;
                double dHMOShare = dSigmaM(nTempHMOidx,nHospStruct,dPrem);
                double dHospShare = dSigmaHM(j,nTempHMOidx,nHospStruct,dPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            
            for (int nTemp=1; nTemp<nH; nTemp++) {                                                  // [N_j'k]
                int nTempHospidx = (j+nTemp)%nH;
                double dHMOShare = dSigmaM(k,nHospStruct,dPrem);
                double dHospShare = dSigmaHM(nTempHospidx,k,nHospStruct,dPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            

            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO k,k',k'', Alternative                   // [N'_jk, N'_jk']
            for (int nTemp=0; nTemp<nM; nTemp++) {
                int nTempHMOidx = (k+nTemp)%nM;
                double dHMOShare = dSigmaM(nTempHMOidx,nHospStruct_New,dNewPrem);
                double dHospShare = dSigmaHM(j,nTempHMOidx,nHospStruct_New,dNewPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }            
            //NUMBER OF PATIENTS SERVED BY HOSPITAL and HMO                                         // [N'_j'k]
            for (int nTemp=1; nTemp<nH; nTemp++) {
                int nTempHospidx = (j+nTemp)%nH;
                double dHMOShare = dSigmaM(k,nHospStruct_New,dNewPrem);
                double dHospShare = dSigmaHM(nTempHospidx,k,nHospStruct_New,dNewPrem);
                double dNMH = nN * dHMOShare * dHospShare;
                oFile << dNMH << "\t";
            }
            
            oFile << dPrem[k] << "\t" << dNewPrem[k] << "\t";                                       // [P_k, P'_k]

            for (int nTemp=0; nTemp<nM; nTemp++) {                                                  // [SigmaM_k, SigmaM_k']
                oFile << dSigmaM((k+nTemp)%nM, nHospStruct, dPrem) << "\t";
            }
            for (int nTemp=0; nTemp<nM; nTemp++) {                                                  // [SigmaM'_k, SigmaM'_k']
                oFile << dSigmaM((k+nTemp)%nM, nHospStruct_New, dNewPrem) << "\t";
            }

            //AllHospCosts                                                                          // [HospCosts: C_j C_j']
            for (int nTemp=0; nTemp<nH; nTemp++) {
                oFile << dCostHospConst[(j+nTemp)%nH] << "\t";
            }            
            
            oFile << dRealizedHospProfits << "\t" << dNewHospProfits << "\t";                       // [RealHospProfits, NewHospProfits]
            oFile << dRealizedHospCosts << "\t" << dNewHospCosts << "\t";                           // [RealHospCosts,NewHospCosts]

            oFile << dPHMO_HospStruct(k,nHospStruct,dTrans) << "\t" << dPHMO_HospStruct(k,nHospStruct_New,dTrans) << "\t"; // [HMOProfits NewHMOProfits]
        } // for nC
        oFile << endl;
    } // for k
    oFile.close();
}



void vCheckIsNEPrems(double dTrans[nMMax][nHMax]) {
    bIsNEPrems = true;
    ofstream oFile("output.txt",ios::app);

    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);    

    double dPrem[nMMax];
    vSetPremiums(dPrem, nHospStruct,dTrans);
    cout << "FOC" << dPremiumFOC(0,dPrem,nHospStruct,dTrans) << " " << dPremiumFOC(1,dPrem,nHospStruct,dTrans) << endl;
    
    double dLambda[nMMax][nRMax];
    double dDelta[nRMax];
    double dPartialSigmaM[nMMax];
    double dHMOTransfers[nMMax];
    double dEqHMOCosts[nMMax];

    for (int r=0; r<nR; r++) {
        for (int k=0; k < nM; k++) {
            double dTemp = 0;
            for (int j=0; j< nH; j++) {
                if (nHospStruct[j]&(int)pow(2.,k)) {
                    dTemp = dTemp + dR_ExpCharHosp[r][j];
                }
            }
            if (dTemp != 0) {
                dLambda[k][r] = exp(dR_CharHMO[r][k]*dCharHMO[k]-(dAlphaPrice*dPrem[k]) + ( dPercentUseHosp*log(dTemp)) );
            } else {
                dLambda[k][r] = exp(dR_CharHMO[r][k]*dCharHMO[k]-(dAlphaPrice*dPrem[k]));
            }
        }
        dDelta[r] = 1;
        for (int k=0; k < nM; k++) {
            dDelta[r] = dDelta[r] + dLambda[k][r];
        }
    }

    cout << "HMO Eq. Costs: \n";
    oFile << "HMO Eq. Costs: \n";

    for (int k=0; k<nM; k++) {
        dHMOTransfers[k]=0;
        for (int j=0; j<nH; j++) {
            dHMOTransfers[k] = dHMOTransfers[k] + dSigmaHM(j,k,nHospStruct,dPrem)* dTrans[k][j];
        }
        dPartialSigmaM[k] = 0;
        for (int r=0; r<nR; r++) {
            dPartialSigmaM[k] = dPartialSigmaM[k] + dR_Shares[r]*dAlphaPrice*(dLambda[k][r]*(dLambda[k][r]-dDelta[r]))/(dDelta[r]*dDelta[r]);
        }
        if (dPartialSigmaM[k] != 0) {
            dEqHMOCosts[k] = dPrem[k] - dHMOTransfers[k] + dSigmaM(k,nHospStruct,dPrem)/dPartialSigmaM[k];
        } else {
            dEqHMOCosts[k] = 0;
        }
        if ((dEqHMOCosts[k] > dPrem[k])||(dEqHMOCosts[k]<0)) bIsNEPrems=false;

        cout << "\t\t" << dEqHMOCosts[k] << "\t \t";
        oFile << "\t\t" << dEqHMOCosts[k] << "\t \t";
        cout << dAlphaPrice*dPrem[k] << " - " << dHMOTransfers[k] << " + " << dSigmaM(k,nHospStruct,dPrem)/dPartialSigmaM[k] << endl;
        oFile << dAlphaPrice*dPrem[k] << " - " << dHMOTransfers[k] << " + " << dSigmaM(k,nHospStruct,dPrem)/dPartialSigmaM[k] << endl;
    }
    cout << endl;
    oFile << endl;
    if (!bIsNEPrems) {    
        cout << "ERROR: No Costs Support Premiums" << endl;
        oFile << "ERROR: No Costs Support Premiums" << endl;
    }
    oFile.close();
}

void vCopyStruct(int nStructCopy[nHMax], int nHospStruct[nHMax]) {
    for (int i=0;i<nHMax;i++) nStructCopy[i] = nHospStruct[i];
}

double dSigmaM (int nMidx, int nHospStruct[nHMax], double dPrem[nMMax]) {
    double dTemp = 0;
    for (int r=0; r<nR; r++) {
        dTemp = dTemp + dR_Shares[r]*dSigmaM_R(nMidx,r,nHospStruct,dPrem);
    }
    return dTemp;
}

double dSigmaHM (int nHidx, int nMidx, int nHospStruct[nHMax],double dPrem[nMMax]) {
    double dTemp = 0;
    for (int r=0; r<nR; r++) {
        dTemp = dTemp + dSigmaR_M(nMidx,r,nHospStruct,dPrem)*dSigmaH_R(nHidx,r,nMidx, nHospStruct,dPrem);
    }
    return dPercentUseHosp*dTemp;
}

double dSigmaH_R (int nHidx, int nRidx, int nMidx, int nHospStruct[nHMax],double dPrem[nMMax]) {
// Calculates share of demographic group nRidx in HMO nMidx who choose Hospital nHidx
    if (nHospStruct[nHidx]&(int)pow(2.,nMidx)) {
        //Total Quality of all hospitals in HMO nMidx
        double dHMOQuality =0;
        for (int ncH=0; ncH < nH; ncH++) {
            if (nHospStruct[ncH]&int(pow(2.,nMidx))) {
                dHMOQuality = dHMOQuality + dR_ExpCharHosp[nRidx][ncH];
            }
        }
        return (dR_ExpCharHosp[nRidx][nHidx]/dHMOQuality);
    } else {
        return 0;
    }
}

double dSigmaM_R (int nMidx, int nRidx, int nHospStruct[nHMax],double dPrem[nMMax]) {
// Calculates share of demographic group nRidx who choose HMO nMidx
    double dHMOQuality[nMMax];
    for (int ncM=0; ncM < nM; ncM++) {
        double dTemp = 0;
        for (int ncH=0; ncH < nH; ncH++) {
            if (nHospStruct[ncH]&(int)pow(2.,ncM)) {
                dTemp = dTemp + dR_ExpCharHosp[nRidx][ncH];
            }
        }
        if (dTemp != 0) {
            dHMOQuality[ncM] = exp(dR_CharHMO[nRidx][ncM]*dCharHMO[ncM]-(dAlphaPrice*dPrem[ncM]) + ( dPercentUseHosp*log(dTemp)) );
        } else {
            dHMOQuality[ncM] = exp(dR_CharHMO[nRidx][ncM]*dCharHMO[ncM]-(dAlphaPrice*dPrem[ncM]));        
        }
    }
    double dNumerator = dHMOQuality[nMidx];
    double dDenominator = 1;
    for (int ncM=0; ncM < nM; ncM++) {
        dDenominator = dDenominator + dHMOQuality[ncM];
    }
    return dNumerator / dDenominator;
}

double dSigmaR_M (int nMidx, int nRidx, int nHospStruct[nHMax],double dPrem[nMMax]) {
// Calculates Demographic Makeup of Ppl who Choose HMO nMidx
    double dSigmaRTimesSigmaMR[nRMax];
    double dTemp = 0;
    for (int r=0; r < nR; r++) {
        dSigmaRTimesSigmaMR[r] = dR_Shares[r] * dSigmaM_R(nMidx,r,nHospStruct,dPrem);
    }
    double dNumerator = dSigmaRTimesSigmaMR[nRidx];
    double dDenominator = 0;
    for (int r=0; r < nR; r++) {
        dDenominator = dDenominator + dSigmaRTimesSigmaMR[r];
    }
    if (dDenominator != 0) {
        return dNumerator / dDenominator;
    } else {
        return 0;
    }
}

void vSetPremiums (double dPrem[nMMax], int nHospStruct[nHMax], double dTrans[nMMax][nHMax]) {
    double dPremInit[nMMax];
    double dPremEnd[nMMax];
    for (int i=0;i<nM;i++) {
        dPremInit[i]=dP[i];
    }
    double dTotalTol = .1;
    double dTotalPremError = 100;
    double dTol = 0;
    while (dTotalPremError > dTotalTol) { 
        for (int i=0;i<nM;i++) {
            dPremEnd[i]=dPremInit[i];
        }
        for (int k=0;k<nM; k++) {
            dPremEnd[k]=dCostHMO[k];
            double dFOC = 100;
            while (dFOC > dTol) {
                dPremEnd[k] = dPremEnd[k] + .01;
                dFOC = (dPremiumFOC(k,dPremEnd,nHospStruct,dTrans));
//                cout.precision(3);
//                cout << dFOC << " ";
            }
            dPremEnd[k] = dPremEnd[k] - .01;
//            cout << endl;
        }
        dTotalPremError = 0;
        for (int i=0;i<nM;i++) {
            dTotalPremError = dTotalPremError + abs(dPremEnd[i]-dPremInit[i]);
            dPremInit[i]=dPremEnd[i];
        }
//        cout << dTotalPremError << " ";
    }
    for (int i=0;i<nM;i++) {
        dPrem[i]=dPremEnd[i];
//        cout << dPrem[i] << " ";
    }
}

double dPremiumFOC(int nMidx, double dPrem[nMMax], int nHospStruct[nHMax], double dTrans[nMMax][nHMax]) {
    double dLambda[nMMax][nRMax];
    double dDelta[nRMax];        
    for (int r=0; r<nR; r++) {
        for (int k=0; k < nM; k++) {
            double dTemp = 0;
            for (int j=0; j< nH; j++) {
                if (nHospStruct[j]&(int)pow(2.,k)) {
                    dTemp = dTemp + dR_ExpCharHosp[r][j];
                }
            }
            if (dTemp != 0) {
                dLambda[k][r] = exp(dR_CharHMO[r][k]*dCharHMO[k]-(dAlphaPrice*dPrem[k]) + ( dPercentUseHosp*log(dTemp)) );
            } else {
                dLambda[k][r] = exp(dR_CharHMO[r][k]*dCharHMO[k]-(dAlphaPrice*dPrem[k]));
            }
        }
        dDelta[r] = 1;
        for (int k=0; k < nM; k++) {
            dDelta[r] = dDelta[r] + dLambda[k][r];
        }
    }
    double dMTransfers=0;
    for (int i=0; i<nH; i++) {
        dMTransfers = dMTransfers +  dSigmaHM(i,nMidx,nHospStruct,dPrem) * dTrans[nMidx][i];
    }

    double dPartialSigmaM = 0;
    for (int r=0; r<nR; r++) {
        dPartialSigmaM = dPartialSigmaM + dR_Shares[r]*dAlphaPrice*(dLambda[nMidx][r]*(dLambda[nMidx][r]-dDelta[r]))/(dDelta[r]*dDelta[r]);
    }    
   
    double dTemp = dPartialSigmaM*(dPrem[nMidx]-dCostHMO[nMidx]-dMTransfers) + dSigmaM(nMidx,nHospStruct,dPrem);
    return dTemp;
}


//=====================================================================================
double dPHosp (int nHidx, double dTrans[nMMax][nHMax]) {
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);    
    return dPHosp_Struct(nHidx, nHospStruct, dTrans);
}

double dPHosp_Struct (int nHidx, int nHospStruct[nHMax], double dTrans[nMMax][nHMax]) {
    double dPrem[nMMax];
    vSetPremiums(dPrem,nHospStruct,dTrans);
    double dHRevenues=0;
    double dHPatients=0;
    double dHMargCosts = 0;
    for (int k=0; k<nM; k++) {
        double dTemp = nN * dSigmaM(k,nHospStruct,dPrem) * dSigmaHM(nHidx,k,nHospStruct,dPrem);
        dHPatients = dHPatients + dTemp;
        dHRevenues = dHRevenues + dTemp * dTrans[k][nHidx];
    }
    dHMargCosts = dCostHospConst[nHidx];
    if (dHPatients > dCapHosp[nHidx]) dHMargCosts = dHMargCosts + pow(dHPatients-dCapHosp[nHidx],dCostHospExp[nHidx]);
    return dHRevenues - dHMargCosts*dHPatients;
}

double dCostsHosp(int nHidx, int nHospStruct[nMMax], double dPrem[nMMax]) {
    double dHPatients=0;
    double dHMargCosts = 0;
    for (int k=0; k<nM; k++) {
        double dTemp = nN * dSigmaM(k,nHospStruct,dPrem) * dSigmaHM(nHidx,k,nHospStruct,dPrem);
        dHPatients = dHPatients + dTemp;
    }
    dHMargCosts = dCostHospConst[nHidx];
    if (dHPatients > dCapHosp[nHidx]) dHMargCosts = dHMargCosts + pow(dHPatients-dCapHosp[nHidx],dCostHospExp[nHidx]);
    return dHMargCosts*dHPatients;
}

double dPHMO (int nMidx, int nHMOStruct[nMMax], double dTrans[nMMax][nHMax]) {
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);    
    double dPrem[nMMax];
    vSetPremiums(dPrem, nHospStruct,dTrans);
    double dMNetPremiums=0;
    double dMTransfers=0;
    for (int i=0; i<nH; i++) {
        double dTemp = nN * dSigmaM(nMidx,nHospStruct,dPrem) * dSigmaHM(i,nMidx,nHospStruct,dPrem);
        dMTransfers = dMTransfers + dTemp * dTrans[nMidx][i];
    }
    dMNetPremiums = nN * dSigmaM(nMidx,nHospStruct,dPrem) * (dPrem[nMidx] - dCostHMO[nMidx]);
    return (dMNetPremiums - dMTransfers);
}

double dPHMO_HospStruct (int nMidx, int nHospStruct[nHMax], double dTrans[nMMax][nHMax]) {
    double dPrem[nMMax];
    vSetPremiums(dPrem, nHospStruct,dTrans);
    double dMNetPremiums=0;
    double dMTransfers=0;
    for (int i=0; i<nH; i++) {
        double dTemp = nN * dSigmaM(nMidx,nHospStruct,dPrem) * dSigmaHM(i,nMidx,nHospStruct,dPrem);
        dMTransfers = dMTransfers + dTemp * dTrans[nMidx][i];
    }
    dMNetPremiums = nN * dSigmaM(nMidx,nHospStruct,dPrem) * (dPrem[nMidx]-dCostHMO[nMidx]);
    return (dMNetPremiums - dMTransfers);
}

void vHMO_GenNE(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]) {
    int nCycles = 0;
    bool bIsCycle=false;
    int nStructHistory[nMMax][nMaxCycles];
    int nStructTemp[nMMax];
    for (int k=0;k<nMMax;k++) {
        nHMOStruct[k]=nS-1;
        nStructTemp[k]=nS-1;
    }
    bool bIsHMONE=false;
    while (!bIsHMONE) {
        bIsHMONE=true;
        for (int k=0;k<nM;k++) nStructTemp[k] = nHMOStruct[k];         // Make a copy of Structure
        vHMO_GenBR(nHMOStruct,dTrans);
        for (int k=0;k<nM;k++) {
            if (nHMOStruct[k]!=nStructTemp[k]) bIsHMONE=false;
        }

        if (nCycles==nMaxCycles) {
            bIsCycle = true;
        } else {
            for (int k=0;k<nM;k++) nStructHistory[k][nCycles] = nHMOStruct[k];
            for (int i=0;i<nCycles-1;i++) {
                bIsCycle=true;
                for (int k=0; k<nM; k++) {
                    if (nStructHistory[k][i]!=nHMOStruct[k]) {
                       bIsCycle=false;
                    }
                }
                if (bIsCycle) {break;}
            }
        }
        if (bIsCycle) {                             // If we can't find a NE via iterative BR procedure, scan the entire space
            cout << "FULL!" << endl;
            vHMO_GenNE_fullsearch(nHMOStruct, dTrans);
            break;
        }
        nCycles++;
    }
}

void vHMO_GenNE_fullsearch(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]) {
    int nStruct_NE[nMMax][nMaxCycles];
    int nNumberOfNEs = 0;
    for (int nC = 0; nC < (int)pow((double)nS,nM); nC++) {
        for (int i=0; i<nM; i++) {
            nHMOStruct[i] = (nC % (int)pow((double)nS,nM-i)) / (int)pow((double)nS,nM-1-i);
        }
        if (bHMO_IsNE(nHMOStruct,dTrans)) {
            for (int i=0;i<nH;i++) nStruct_NE[i][nNumberOfNEs] = nHMOStruct[i];
            nNumberOfNEs++;
        }
    }
    if (nNumberOfNEs == 0) {
        bIsHMOCycling = true;
        cout << "WARNING - NO PURE STRATEGIES EXIST FOR THIS PARTICULAR SET OF TRANSFERS: \n";
        vPrintTransfers(dTrans);
    } else {
        // Equilibrium Selection : Choose the one that maximizes HMO Profits
        int nTempIdx=0;
        double dTempProfits=0;
        for (int j=0; j<nNumberOfNEs; j++) {
            double dTempProfits2=0;
            int nStructTemp[nMMax];
            for (int i=0;i<nM;i++) nStructTemp[i] = nStruct_NE[i][j];
            for (int i=0;i<nM;i++) dTempProfits2 = dTempProfits2 + dPHMO(i,nHMOStruct,dTrans);
            if (dTempProfits2>dTempProfits) nTempIdx = j;
        }

        for (int i=0;i<nM;i++) nHMOStruct[i] = nStruct_NE[i][nNumberOfNEs-1];
//        cout << "(" << nTempIdx+1 << ":" << nNumberOfNEs <<") " << endl;
//        ofstream oFile("output.txt",ios::app);
//        oFile << "(" << nTempIdx+1 << ":" << nNumberOfNEs <<") ";
//        oFile.close();         
    }
}

bool bHMO_IsNE(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]) {
    int nStructTemp[nMMax];
    bool bIsHMONE=true;
    for (int i=0;i<nM;i++) nStructTemp[i] = nHMOStruct[i];         // Make a copy of Structure
    vHMO_GenBR(nHMOStruct,dTrans);
    for (int i=0;i<nH;i++) {
        if (nHMOStruct[i]!=nStructTemp[i]) {
            bIsHMONE=false;
        }
    }
    return bIsHMONE;
}

void vHMO_GenBR(int nHMOStruct[nMMax],double dTrans[nMMax][nHMax]) {
    for (int i=0;i<nM;i++) {
        int nStratBR = nHMOStruct[i];
        double dCurrentProfits = dPHMO(i,nHMOStruct,dTrans);
        for (int nStrat=0; nStrat < nS; nStrat++) {
            nHMOStruct[i] = nStrat;
            double dNewProfits = dPHMO(i,nHMOStruct,dTrans);
            if (dNewProfits > dCurrentProfits) {
                nStratBR = nStrat;
                dCurrentProfits = dNewProfits;
            }
        }
        nHMOStruct[i] = nStratBR;
    }
}

void vConvertStruct_HMO2Hosp(int nHMOStruct[nMMax], int nHospStruct[nHMax]) {
    int nIndustryStructure[nMMax][nHMax];
    for (int k=0; k<nM; k++) {
        for (int j=0; j<nH; j++) {
            if (nHMOStruct[k]&(int)pow(2.,j)) {
                nIndustryStructure[k][j]=1;                
            } else {
                nIndustryStructure[k][j]=0;
            }
        }   
    }    

    for (int j=0; j<nH; j++) {
        nHospStruct[j] = 0;
        for (int k=0; k<nM; k++) {
            if (nIndustryStructure[k][j]) nHospStruct[j]=nHospStruct[j] + (int)pow(2.,k);
        }
    }
}

void vChangeStruct_HM(int nHidx, int nMidx, int nHospStruct[nHMax], int nHospStructNew[nHMax]) {
    int nIndustryStructure[nMMax][nHMax];
    for (int j=0; j<nH; j++) {
        for (int k=0; k<nM; k++) {
            if (nHospStruct[j]&(int)pow(2.,k)) {
                nIndustryStructure[k][j]=1;                
            } else {
                nIndustryStructure[k][j]=0;
            }
        }   
    }
    
    nIndustryStructure[nMidx][nHidx] = (nIndustryStructure[nMidx][nHidx] + 1) % 2;

    for (int j=0; j<nH; j++) {
        nHospStructNew[j] = 0;
        for (int k=0; k<nM; k++) {
            if (nIndustryStructure[k][j]) nHospStructNew[j] = nHospStructNew[j] + (int)pow(2.,k);
        }
    }
}

void vChangeStruct_HM2(int nHidx, int nMidx, int nHospStruct[nHMax], int nHospStructNew[nHMax]) {
    int nIndustryStructure[nMMax][nHMax];
    for (int j=0; j<nH; j++) {
        for (int k=0; k<nM; k++) {
            if (nHospStruct[j]&(int)pow(2.,k)) {
                nIndustryStructure[k][j]=1;                
            } else {
                nIndustryStructure[k][j]=0;
            }
        }   
    }
    for (int j=0; j<nH; j++) {
        nIndustryStructure[nMidx][j] = (nIndustryStructure[nMidx][j] + 1) % 2;
    }
    for (int j=0; j<nH; j++) {
        nHospStructNew[j] = 0;
        for (int k=0; k<nM; k++) {
            if (nIndustryStructure[k][j]) nHospStructNew[j] = nHospStructNew[j] + (int)pow(2.,k);
        }
    }
}



void vCopyTransfers(double dTransCopy[nMMax][nHMax], double dTrans[nMMax][nHMax]) {
    for (int k=0;k<nM;k++){
        for (int i=0;i<nH;i++){
            dTransCopy[k][i] = dTrans[k][i];
        }
    }
}

void vPrintTransfers(double dTrans[nMMax][nHMax]) {
    ofstream oFile("output.txt",ios::app);    
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout.precision(1);
    cout << "  ";
    
    for (int j=0;j<nH;j++){
        cout << "(";
        for (int k=0;k<nM;k++){
            cout << dTrans[k][j];
            if (k!=(nM-1)) {
                cout <<"  ";
            }
        }
        cout << ")";
        if (j!=(nH-1)) {
            cout <<" , ";
        }
    }
    cout << endl;
}

void vPrintTransfers_File(double dTrans[nMMax][nHMax]) {
    ofstream oFile("output.txt",ios::app);
    oFile.setf(ios::fixed);
    oFile.setf(ios::showpoint);
    oFile.precision(3);
    oFile << "  ";
    
    for (int j=0;j<nH;j++){
        oFile << "(";
        for (int k=0;k<nM;k++){
            oFile << dTrans[k][j];
            if (k!=(nM-1)) {
                oFile <<"  ";
            }
        }
        oFile << ")";
        if (j!=(nH-1)) {
              oFile <<" , ";
        }
    }
    oFile << endl;    
    oFile.close();
}

void vPrintHMOStructure_FromTrans(double dTrans[nMMax][nHMax]) {
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    vPrintHMOStructure(nHMOStruct);
}

void vPrintHMOStructure(int nHMOStruct[nHMax]) {
    cout << "\t\t\t\t\t\t[ ";
    for (int k=0;k<nM;k++) cout << nHMOStruct[k] << " ";
    cout << "]" << endl;    

    ofstream oFile("output.txt",ios::app);    
    oFile << "\t\t\t\t\t\t[ ";
    for (int k=0;k<nM;k++) oFile << nHMOStruct[k] << " ";
    oFile << "]" << endl;    
}

void vPrintHospStructure_FromTrans(double dTrans[nMMax][nHMax]) {
    int nHMOStruct[nMMax];
    vHMO_GenNE(nHMOStruct,dTrans);
    int nHospStruct[nHMax];
    vConvertStruct_HMO2Hosp(nHMOStruct,nHospStruct);        
    vPrintHospStructure(nHospStruct);
}

void vPrintHospStructure(int nHospStruct[nHMax]) {
    cout << "\t\t\t\t\t\t[ ";
    for (int j=0;j<nH;j++) cout << nHospStruct[j] << " ";
    cout << "]" << endl;    

    ofstream oFile("output.txt",ios::app);    
    oFile << "\t\t\t\t\t\t[ ";
    for (int j=0;j<nH;j++) oFile << nHospStruct[j] << " ";
    oFile << "]" << endl;    
}

int nIndustryStructure_FromHMOStruct(int nHMOStruct[nMMax]) {
    int nIndustryStruct = 0;
    for (int k=0;k<nM;k++) {    
        nIndustryStruct = nIndustryStruct + nHMOStruct[k]*(int)pow( (double)nS,nM-1-k );
    }
    return nIndustryStruct;
}

void vPrintIndustryStructure_FromHMOStruct(int nHMOStruct[nMMax]) {
    int nIndustryStructure[nMMax][nHMax];
    for (int k=0; k<nM; k++) {
        for (int j=0; j<nH; j++) {
            if (nHMOStruct[k]&(int)pow(2.,j)) {
                nIndustryStructure[k][j]=1;                
            } else {
                nIndustryStructure[k][j]=0;
            }
        }   
    }
    cout << "Industry Structure:\n";
    for (int k=0; k<nM; k++) {
        cout << "\t HMO " << k+1 << ":\t";
        for (int j=0; j<nH; j++) {
            cout << nIndustryStructure[k][j] << " ";
        }
        cout << endl;
    }
}

void vPrintIndustryStructure_FromHMOStruct_File(int nHMOStruct[nMMax]) {
    ofstream oFile("output.txt",ios::app);
    int nIndustryStructure[nMMax][nHMax];
    for (int k=0; k<nM; k++) {
        for (int j=0; j<nH; j++) {
            if (nHMOStruct[k]&(int)pow(2.,j)) {
                nIndustryStructure[k][j]=1;                
            } else {
                nIndustryStructure[k][j]=0;
            }
        }   
    }
    oFile << "Industry Structure:\n";
    for (int k=0; k<nM; k++) {
        oFile << "\t HMO " << k+1 << ":\t";
        for (int j=0; j<nH; j++) {
            oFile << nIndustryStructure[k][j] << " ";
        }
        oFile << endl;
    }
    oFile.close();
}
 

//******************************************************************************
double normal_01_cdf ( double x )

//******************************************************************************
//
//  Purpose:
//
//    NORMAL_01_CDF evaluates the Normal 01 CDF.
//
//  Modified:
//
//    10 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference: 
//
//    A G Adams,
//    Areas Under the Normal Curve,
//    Algorithm 39, 
//    Computer j., 
//    Volume 12, pages 197-198, 1969.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
  double a1 = 0.398942280444E+00;
  double a2 = 0.399903438504E+00;
  double a3 = 5.75885480458E+00;
  double a4 = 29.8213557808E+00;
  double a5 = 2.62433121679E+00;
  double a6 = 48.6959930692E+00;
  double a7 = 5.92885724438E+00;
  double b0 = 0.398942280385E+00;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302E+00;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364E+00;
  double b5 = 0.151679116635E+00;
  double b6 = 5.29330324926E+00;
  double b7 = 4.8385912808E+00;
  double b8 = 15.1508972451E+00;
  double b9 = 0.742380924027E+00;
  double b10 = 30.789933034E+00;
  double b11 = 3.99019417011E+00;
  double cdf;
  double q;
  double y;
//
//  |X| <= 1.28.
//
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 
      + a6 / ( y + a7 ) ) ) );
//
//  1.28 < |X| <= 12.7
//
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1 
      + b2 / ( fabs ( x ) + b3 
      + b4 / ( fabs ( x ) - b5 
      + b6 / ( fabs ( x ) + b7 
      - b8 / ( fabs ( x ) + b9 
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
//
//  12.7 < |X|
//
  }
  else
  {
    q = 0.0;
  }
//
//  Take account of negative X.
//
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}
//**********************************************************************

double normal_01_cdf_inv ( double p )

//**********************************************************************
//
//  Purpose:
//
//    NORMAL_01_CDF_INV inverts the standard normal CDF.
//
//  Discussion:
//
//    The result is accurate to about 1 part in 10**16.
//
//  Modified:
//
//    27 December 2004
//
//  Author:
//
//    Michael Wichura.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Michael Wichura,
//    The Percentage Points of the Normal Distribution,
//    Algorithm AS 241,
//    Applied Statistics,
//    Volume 37, Number 3, pages 477-484, 1988.
//
//  Parameters:
//
//    Input, real ( kind = 8 ) P, the value of the cumulative probability 
//    densitity function.  0 < P < 1.  If P is outside this range, an
//    "infinite" value is returned.
//
//    Output, real ( kind = 8 ) NORMAL_01_CDF_INVERSE, the normal deviate value 
//    with the property that the probability of a standard normal deviate being 
//    less than or equal to this value is P.
//
{
  double a[8] = {
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
  double b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
  double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770e-1,
    2.27238449892691845833e-2,  7.74545014278341407640e-4 };
  double const1 = 0.180625;
  double const2 = 1.6;
  double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550e-1,
    1.48103976427480074590e-1,  1.51986665636164571966e-2,
    5.47593808499534494600e-4,  1.05075007164441684324e-9 };
  double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230e-1,
    2.65321895265761230930e-2,  1.24266094738807843860e-3,
    2.71155556874348757815e-5,  2.01033439929228813265e-7 };
  double f[8] = {
    1.0,                        5.99832206555887937690e-1,
    1.36929880922735805310e-1,  1.48753612908506148525e-2,
    7.86869131145613259100e-4,  1.84631831751005468180e-5,
    1.42151175831644588870e-7,  2.04426310338993978564e-15 };
  double q;
  double r;
  double split1 = 0.425;
  double split2 = 5.0;
  double value;

  if ( p <= 0.0 )
  {
    value = -d_huge ( );
    return value;
  }

  if ( 1.0 <= p )
  {
    value = d_huge ( );
    return value;
  }

  q = p - 0.5;

  if ( fabs ( q ) <= split1 )
  {
    r = const1 - q * q;
    value = q * dpoly_value ( 8, a, r ) / dpoly_value ( 8, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = -1.0;
      exit ( 1 );
    }

    r = sqrt ( -log ( r ) );

    if ( r <= split2 )
    {
      r = r - const2;
      value = dpoly_value ( 8, c, r ) / dpoly_value ( 8, d, r );
     }
     else
     {
       r = r - split2;
       value = dpoly_value ( 8, e, r ) / dpoly_value ( 8, f, r );
    }

    if ( q < 0.0 )
    {
      value = -value;
    }

  }

  return value;
}

//**********************************************************************

double d_huge ( void )

//******************************************************************************
//
//  Purpose:
//
//    D_HUGE returns a "huge" real value, usually the largest legal real.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_HUGE, a "huge" real value.
//
{
  return HUGE_VAL;
}

//**********************************************************************

double dpoly_value ( int n, double a[], double x )

//**********************************************************************
//
//  Purpose:
//
//    DPOLY_VALUE evaluates a double precision polynomial.
//
//  Discussion:
//
//    For sanity's sake, the value of N indicates the NUMBER of 
//    coefficients, or more precisely, the ORDER of the polynomial,
//    rather than the DEGREE of the polynomial.  The two quantities
//    differ by 1, but cause a great deal of confusion.
//
//    Given N and A, the form of the polynomial is:
//
//      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Input, double A[N], the coefficients of the polynomial.
//    A[0] is the constant term.
//
//    Input, double X, the point at which the polynomial is to be evaluated.
//
//    Output, double DPOLY_VALUE, the value of the polynomial at X.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}

//**********************************************************************

double normal_01_sample ( int *seed )

//**********************************************************************
//
//  Purpose:
//
//    NORMAL_01_SAMPLE samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has 
//    mean 0 and standard deviation 1.
//
//  Method:
//
//    The Box-Muller method is used, which is efficient, but 
//    generates two values at a time.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double NORMAL_01_SAMPLE, a normally distributed random value.
//
{
 # define PI 3.141592653589793

  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = d_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = d_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * PI * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * PI * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
 # undef PI
}

//******************************************************************************

int get_seed ( void )

//******************************************************************************
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
  time_t clock;
  int i;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  clock = time ( &tloc );
  lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
//
//  Move Hours to 0, 1, ..., 11
//
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
  seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,IMAX].
//
  seed = ( int ) 
    ( ( ( double ) seed )
    * ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
#undef I_MAX
}

//******************************************************************************

double d_uniform_01 ( int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D_UNIFORM_01 is a portable pseudorandom number generator.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double D_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
