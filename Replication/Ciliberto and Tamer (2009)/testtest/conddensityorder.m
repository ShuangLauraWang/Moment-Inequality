function prob=conddensityorder()

filename='marketdata.csv';
marketdata=dataset ('File', filename, 'Delimiter',',');
datairline=marketdata;

%MarketData.Properties.VarNames

global k Y rowX total


%To simplify the model, I will first try to include only two firms, Delta
%&AA, and take the entry decisions of other firms as given. Also, I will
%only include airline presence and cost as explanory variables. 
%(maybe include market size later)

% To discretize the Xs, I first plot their densities 
% [f1_DL,xi1_DL] = ksdensity(MarketData.marketpresenceDL);       
% [f1_AA,xi1_AA] = ksdensity(MarketData.marketpresenceAA);
% figure
% plot(xi1_DL,f1_DL,xi1_AA,f1_AA); %quite normally distributed
%                                  quartiles
%                                  %DL 0.41, 0.55, 0.68
%                                  %AA 0.29, 0.42, 0.55

% Since the marketpresence is "share", it is distributed between 0 to 1. I
% categorize them into 0.2 intervals.
% 
% [f2_DL,xi2_DL] = ksdensity(MarketData.mindistancefromhubDL);       
% [f2_AA,xi2_AA] = ksdensity(MarketData.mindistancefromhubAA);
% figure
% plot(xi2_DL,f2_DL,xi2_AA,f2_AA); %skewed
%                                  %quartiles
%                                  %DL 0.005, 0.051, 0.269
%                                  %AA 0.016, 0.170, 0.812 (4*4)^2=256

%a simple frequency estimator 
y=PossibleOutcome(2,2)-1;

%IY is a matrix, the ith column indicates wether the market is of ith
%market structure
IY=zeros(rowX,2^k);

for i=1:size(y,1)
IY(:,i)=ismember(Y,y(i,:),'rows');
end

%Find out the market structure of each market
[row,col]=find(IY);
YCat=[row col];
YCat=sortrows(YCat,1);

clear row col

%convert IY into a numeric matrix
IY=double(IY);


p = 0:.25:1;
% breaks = quantile(datairline.marketdistance,p);
% datairline.mdQ = ordinal(datairline.marketdistance,{'mdQ1','mdQ2',...
%     'mdQ3','mdQ4'},[],breaks);
% 
% breaks = quantile(datairline.fromcenterdistance,p);
% datairline.cdQ = ordinal(datairline.fromcenterdistance,{'cdQ1','cdQ2',...
%     'cdQ3','cdQ4'},[],breaks);
% 
% breaks = quantile(datairline.mindistance,p);
% datairline.caQ = ordinal(datairline.mindistance,{'caQ1','caQ2',...
%     'caQ3','caQ4'},[],breaks);
% 
% breaks = quantile(datairline.changeincmarket,p);
% datairline.igQ = ordinal(datairline.changeincmarket,{'igQ1','igQ2',...
%     'igQ3','igQ4'},[],breaks);
% 
% breaks = quantile(datairline.percapitaincmarket,p);
% datairline.piQ = ordinal(datairline.percapitaincmarket,{'piQ1','piQ2',...
%     'piQ3','piQ4'},[],breaks);
% 
% breaks = quantile(datairline.marketsize,p);
% datairline.msQ = ordinal(datairline.marketsize,{'msQ1','msQ2',...
%     'msQ3','msQ4'},[],breaks);

breaks = quantile(datairline.marketpresenceAA,p);
datairline.mpAAQ = ordinal(datairline.marketpresenceAA,{'mpAAQ1','mpAAQ2',...
    'mpAAQ3','mpAAQ4'},[],breaks);
               
breaks = quantile(datairline.marketpresenceDL,p);    
datairline.mpDLQ = ordinal(datairline.marketpresenceDL,{'mpDLQ1','mpDLQ2',...
    'mpDLQ3','mpDLQ4'},[],breaks);    

breaks = quantile(datairline.mindistancefromhubAA,p);    
datairline.mdfhAAQ = ordinal(datairline.mindistancefromhubAA,{'mdfhAAQ1',...
    'mdfhAAQ2','mdfhAAQ3','mdfhAAQ4'},[],breaks);  


breaks = quantile(datairline.mindistancefromhubDL,p);    
datairline.mdfhDLQ = ordinal(datairline.mindistancefromhubDL,{'mdfhDLQ1',...
    'mdfhDLQ2','mdfhDLQ3','mdfhDLQ4'},[],breaks);   

X=double(datairline(:,{'mpAAQ','mpDLQ', 'mdfhAAQ','mdfhDLQ'}));
x=PossibleOutcome(4,4);

for i=1:size(x,1)
IX(:,i)=ismember(X,x(i,:),'rows');
end

[row,col]=find(IX);
XCat=[row col];
XCat=sortrows(XCat,1);

IX=double(IX);

xy=IY'*IX;
sumx=sum(IX,1);
PyX=xy./repmat(sumx,2^k,1); %the probability of y(a certain type of market
                           %structure) when X falls into a certian category

%dimension of prob, 2742*64??????????                           
prob=zeros(rowX,total);

for i=1:rowX
    prob(i,:)=PyX(:,XCat(i,2))';
end

end

