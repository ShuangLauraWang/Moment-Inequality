clear;
warning off all;
randn('state',0);

load XALL;

size(x1ALL);
size(x2ALL);
TotalObs = size(x1ALL,1);
TotalMkts = TotalObs / 4;
nMonteDraws = 600;
SegmentSize = 5540;
load Estimates_TrueBetaNu2s_072709.mat;

RandomDraws = zeros(SegmentSize/4,nMonteDraws);
ChunkSize=10;
for iMonteDraw = 1:nMonteDraws
    for i=1:SegmentSize/4
        RandomDraws(i,iMonteDraw) = floor(rand*TotalMkts) + 1;
        ErrFlag = 1;
        while ErrFlag
            ErrFlag = 0;
            if sum(RandomDraws(i,iMonteDraw)==RandomDraws(1:i-1,iMonteDraw))>0
                ErrFlag = 1;
                RandomDraws(i,iMonteDraw) = floor(rand*TotalMkts) + 1;                
            end;
        end;        
    end;    
end;

SDMultipleCost = .25^.5;
SDMultiplePop = .05^.5;
SD1 = 9;
SD2 = 300;
SD_HospMeasurementErrors = SDMultipleCost*SD1;
SD_NBarMeasurementErrors = SDMultiplePop*SD2;

HospMeasurementErrorsAll = normrnd(0,SD_HospMeasurementErrors, TotalObs/2, 1);
NBarMeasurementErrorsAll = normrnd(0,SD_NBarMeasurementErrors, TotalObs/4, 1);

for Chunk=(1):(nMonteDraws/ChunkSize)
    clear XALL;
    disp(Chunk);
    for i = 1:ChunkSize
        %disp(i);
        XALL(i).x1 = zeros(SegmentSize, size(x1ALL,2));
        XALL(i).x2 = zeros(SegmentSize, size(x2ALL,2));   
        XALL(i).HospMeas = zeros(SegmentSize/2, 1);
        XALL(i).NBarMeas = zeros(SegmentSize/4, 1);  
        XALL(i).RealNu2s_MktDraws = zeros(SegmentSize/4, 4);
        for j = 1:SegmentSize/4
            XALL(i).x1(4*(j-1)+1:4*(j-1)+4,:) = x1ALL(4*(RandomDraws(j,(i-1)*Chunk+i)-1)+1:4*(RandomDraws(j,(i-1)*Chunk+i)-1)+4,:);
            XALL(i).x2(4*(j-1)+1:4*(j-1)+4,:) = x2ALL(4*(RandomDraws(j,(i-1)*Chunk+i)-1)+1:4*(RandomDraws(j,(i-1)*Chunk+i)-1)+4,:);        
            XALL(i).HospMeas( 2*(j-1)+1:2*(j-1)+2, 1) = HospMeasurementErrorsAll(2*(RandomDraws(j,(i-1)*Chunk+i)-1)+1:2*(RandomDraws(j,(i-1)*Chunk+i)-1)+2,1);
            XALL(i).NBarMeas(j,1) = NBarMeasurementErrorsAll(RandomDraws(j,(i-1)*Chunk+i),1);
            XALL(i).RealNu2s_MktDraws(j,:) = RealNu2s(RandomDraws(j,(i-1)*Chunk+i),:);
        end;    
    end;

    code=['save DataChunk_072709_' num2str(Chunk) ' XALL HospMeasurementErrorsAll NBarMeasurementErrorsAll'];
    eval(code);
end;