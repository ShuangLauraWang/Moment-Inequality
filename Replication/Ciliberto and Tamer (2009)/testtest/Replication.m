
filename='marketdata.csv';
MarketData=dataset ('File', filename, 'Delimiter',',');

MarketData.Properties.VarNames

global k n
k=6;
n=size(MarketData,2);

%To simplify the model, I will first try to include only two firms, Delta
%&AA, and take the entry decisions of other firms as given. Also, I will
%only include airline presence and cost as explanory variables. 
%(maybe include market size later)

% To discretize the Xs, I first plot their densities 
[f1_DL,xi1_DL] = ksdensity(MarketData.marketpresenceDL);       
[f1_AA,xi1_AA] = ksdensity(MarketData.marketpresenceAA);
figure
plot(xi1_DL,f1_DL,xi1_AA,f1_AA); %quite normally distributed
                                 %quartiles
                                 %DL 0.41, 0.55, 0.68
                                 %AA 0.29, 0.42, 0.55

% Since the marketpresence is "share", it is distributed between 0 to 1. I
% categorize them into 0.2 intervals.

[f2_DL,xi2_DL] = ksdensity(MarketData.mindistancefromhubDL);       
[f2_AA,xi2_AA] = ksdensity(MarketData.mindistancefromhubAA);
figure
plot(xi2_DL,f2_DL,xi2_AA,f2_AA); %skewed
                                 %quartiles
                                 %DL 0.005, 0.051, 0.269
                                 %AA 0.016, 0.170, 0.812 (4*4)^2=256

%a simple frequency estimator 
entry=double(MarketData(:,2:7));
y=PossibleOutcome(2,6)-1;

%IY is a matrix 
for i=1:size(y,1)
IY(:,i)=ismember(entry,y(i,:),'rows');
end

IY=double(IY);

covariates=double([MarketData.marketpresenceAA MarketData.marketpresenceDL...
    MarketData.mindistancefromhubAA MarketData.mindistancefromhubDL]);



MarketData.Properties.VarNames


p = 0:.25:1;
breaks = quantile(MarketData.marketpresenceAA,p);
MarketData.mpAAQ = ordinal(MarketData.marketpresenceAA,{'mpAAQ1','mpAAQ2',...
    'mpAAQ3','mpAAQ4'},[],breaks);
               

breaks = quantile(MarketData.marketpresenceDL,p);    
MarketData.mpDLQ = ordinal(MarketData.marketpresenceDL,{'mpDLQ1','mpDLQ2',...
    'mpDLQ3','mpDLQ4'},[],breaks);    

breaks = quantile(MarketData.mindistancefromhubAA,p);    
MarketData.mdfhAAQ = ordinal(MarketData.mindistancefromhubAA,{'mdfhAAQ1',...
    'mdfhAAQ2','mdfhAAQ3','mdfhAAQ4'},[],breaks);  


breaks = quantile(MarketData.mindistancefromhubDL,p);    
MarketData.mdfhDLQ = ordinal(MarketData.mindistancefromhubDL,{'mdfhDLQ1',...
    'mdfhDLQ2','mdfhDLQ3','mdfhDLQ4'},[],breaks);   

XCat=double(MarketData(:,{'mpAAQ','mpDLQ', 'mdfhAAQ','mdfhDLQ'}));
x=PossibleOutcome(4,4);

for i=1:size(x,1)
IX(:,i)=ismember(XCat,x(i,:),'rows');
end

IX=double(IX);

xy=IY'*IX;
sumx=sum(IX,1);
PyX=xy./repmat(sumx,64,1);

%--------------------------------------------------------------------------%
% (4) SET THE NUMBER OF SIMULATIONS                                        %
%--------------------------------------------------------------------------%
                                                                           %
global r                                                                   %
r=100;

%--------------------------------------------------------------------------%
% (7) DEFINE A COUNTING VARIABLE TO KEEP TRACK OF HOW MANY ITERATIONS HAVE %
% BEEN EXECUTED                                                            %
%--------------------------------------------------------------------------%
                                                                           %
global iteration                                                           %
iteration=0;    



%--------------------------------------------------------------------------%
% (9) CONSTRUCTS INDEXES THAT PERMIT TO AVOID HAVING TO LOOP OVER MARKETS  %
%--------------------------------------------------------------------------%
%

X=covariates;
index=x;

global total
total=size(index,1);

global indexheter                                                          %
tempheter=reshape(1:k*size(X,1),size(X,1),k);                              %
indexheter=zeros(size(X,1)*total,k);                                       %
for j=1:size(X,1)                                                          %
    indexheter((j-1)*total+1:(j-1)*total+total,:)=repmat(tempheter(j,:),total,1);                                                       %
end                                                                        %
clear tempheter                                                            %
                                                                           %
global repindex                                                            %
repindex=repmat(index,size(X,1),1);                                        %
                                                                           %
global diffdummies                                                         %
diffdummies=repmat(reshape(1:total*k,total,k),size(X,1),1);                %

%--------------------------------------------------------------------------%
% (10) CREATE A MATLAB FILE THAT STORES THE RESULTS FROM EACH ITERATION    %
%--------------------------------------------------------------------------%
                                                                           %
oldresultsasa=[0,0,iteration,r,param0'];                                   %
save(['oldresultsasa',sim_num],'oldresultsasa')                            %
                                                  
