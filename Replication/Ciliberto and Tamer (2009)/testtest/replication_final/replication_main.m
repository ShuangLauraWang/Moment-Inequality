%Edited by Shuang Wang for EC711 Final Project
%March 18, 2017

clear;
clc;


%--------------------------------------------------------------------------%
% (1) READS IN THE INITIAL VALUES OF THE PARAMETERS AND AN INDICATOR THAT  %
% CAN BE ASSOCIATED TO THEM                                                %
%--------------------------------------------------------------------------%
                                                                           %
start_values=textread('start_values');                                     % This is a matrix, whose first row is an indicator that is associated
                                                                           % with a particular starting value for the parameters and the second
                                                                           % row is the parameter vector for the starting values.
                                                                           % Dimension: 2 X (# parameters)
                                                                           %
global sim_num                                                             %
sim_num=num2str(start_values(1,1));                                        % Records particular starting value.
                                                                           %
param0 = [start_values(2,1) start_values(2,10:17)]';


%--------------------------------------------------------------------------%
% (2) LOADS DATA                                                           %
%--------------------------------------------------------------------------%
filename='marketdata.csv';
marketdata=dataset ('File', filename, 'Delimiter',',');

datairline=marketdata;

%--------------------------------------------------------------------------%
% (3) DEFINE VARIABLE COUNTING THE NUMBER OF FIRMS, THE MATRIX OF POSSIBLE %
% MARKET STRUCTURES                                                        %
%--------------------------------------------------------------------------%
                                                                           %
global index total k                                                       %
k=2;                                                                       % Number of firms. [Here I only take the entry decisions of AA
                                                                           % and Delta as exogeneous]
                                                                           %
index=makeindex(2);                                                        % Matrix of possible market structures.
                                                                           % Dimension: (2^#firms) X #firms
                                                                           %
total=size(index,1);                                                       % This is the number of possible market structures = 2^#firms.

%--------------------------------------------------------------------------%
% (4) SET THE NUMBER OF SIMULATIONS                                        %
%--------------------------------------------------------------------------%
                                                                           %
global r                                                                   %
r=100;           

%--------------------------------------------------------------------------%
% (5) DEFINE VARIABLES ASSOCIATED WITH DATA                                %
%--------------------------------------------------------------------------%


global X Y rowX colX otherentry                                            %

Xheter=double(datairline(:,{'marketpresenceAA','marketpresenceDL'...
    'mindistancefromhubAA', 'mindistancefromhubDL'}));                     % Exogenous regressors.

otherentry=double(datairline(:,{'airlineUA','airlineAL','airlineLCC',...
    'airlineWN'}));                                                        % [The entry decisions of UA, MA, LCC & WN 
                                                                           % are taken as given].

X=[Xheter otherentry];                                                     % Dimension: # markets X # regressors
                                                                           %
colX=size(X,2);                                                            % Number of exogenous regressors.
                                                                           %
Y=double(datairline(:,{'airlineAA','airlineDL'})); 

rowX=size(X,1); 

%--------------------------------------------------------------------------%
% (6) FIRST STAGE EMPIRICAL PROBABILITIES                                  %
%--------------------------------------------------------------------------%
                                                                           %
global prob                                                                %                                              % Loads empirical probabilities.
prob=conddensityorder;  


%--------------------------------------------------------------------------%
% (7) DEFINE A COUNTING VARIABLE TO KEEP TRACK OF HOW MANY ITERATIONS HAVE %
% BEEN EXECUTED                                                            %
%--------------------------------------------------------------------------%
                                                                           %
global iteration                                                           %
iteration=0;                                                               %
                                                                           %
%--------------------------------------------------------------------------%
% (8) LOAD THE SIMULATED UNOBSERVABLES                                     %
%--------------------------------------------------------------------------%
                                                                           %
global epsi                                                                %
%load firmmarketerrs.raw                                                   %
                                                                           %
%epsi=firmmarketerrs(:,1:r*k);                                             % Firm specific unobservables. We draw outside of this program (in Stata).
rng(1);
epsi=normrnd(0,1,rowX,r*k);                                                % # markets X (# simulations * #firms)
  

%epsimarket=firmmarketerrs(:,r*k+1:r*k+r);                                 % Market specific unobservables.
rng(2);
epsimarket=normrnd(0,1,rowX,r);                                            % Dimension: # markets X # simulations
                                                                           %
epsimarket=kron(epsimarket,ones(1,k));                                     % Assigns the market unobservable to each of the k firms.
                                                                           % Dimension: # markets X (# simulations * #firms)
                                                                           %
epsi=epsi+epsimarket;                                                      % Adds market and firm unobservables.
                                                                           %
                                                                          
                                                                       
market=char(datairline.market);                                            % [Generate airport specific errors, one for each end]
origin=market(:,1:3);
dest=market(:,end-2:end);

airport=unique([origin;dest],'rows');
nrow=size(airport,1);

rng(3)
AirRndErr=normrnd(0,1,nrow,r);

epsiairport=zeros(rowX,r);

for i=1:nrow
    for j=1:rowX
        if origin(j,:)==airport(i,:)
            epsiairport(j,:)=epsiairport(j,:)+AirRndErr(i,:);
        end
        if dest(j,:)==airport(i,:)
            epsiairport(j,:)=epsiairport(j,:)+AirRndErr(i,:);
        end
    end
    
end

epsiairport=kron(epsiairport,ones(1,k));                                   %[Dimension: # markets X (# simulations * #firms)]

epsi=epsi+epsiairport;                                                     %[Sum all unobservables up]

clear nrow

%--------------------------------------------------------------------------%
% (9) CONSTRUCTS INDEXES THAT PERMIT TO AVOID HAVING TO LOOP OVER MARKETS  %
%--------------------------------------------------------------------------%
                                                                           %
global indexheter                                                          %
tempheter=reshape(1:k*size(X,1),size(X,1),k);                              %
indexheter=zeros(size(X,1)*total,k);                                       %
for j=1:size(X,1)                                                          %
    indexheter((j-1)*total+1:(j-1)*total+total,:)=repmat(tempheter(j,:),total,1);                                                       %
end                                                                        %
clear tempheter                                                            %
                                                                           %
global repindex 
repindex=repmat(index,size(X,1),1);                                        %
                                                                           %
global diffdummies                                                         %
diffdummies=repmat(reshape(1:total*k,total,k),size(X,1),1);                %


                                                                           %
%--------------------------------------------------------------------------%
% (10) CREATE A MATLAB FILE THAT STORES THE RESULTS FROM EACH ITERATION    %
%--------------------------------------------------------------------------%
                                                                           %
oldresultsasa=[0,0,iteration,r,param0'];                                   %
save(['oldresultsasa',sim_num],'oldresultsasa') 

%--------------------------------------------------------------------------%
% (11) SIMULATED ANNEALING                                                 %
%--------------------------------------------------------------------------%
                                                                           %
lb = ones(1,9)*(-30);                                                      % Minimum values that the parameters can take.
ub = ones(1,9)*30;                                                         % Maximum values that the parameters can take.
                                                                           %
[param,fval] = simulannealbnd(@replication_simuhete,param0,lb, ub);        % Here, it is useful to try different temperatures, different schedules
                                                                           % for the changes in temperatures, and so on, especially when starting 
                                                                           % the miminization process.

%--------------------------------------------------------------------------%
% (12) CONFIDENCE REGION                                            %
%--------------------------------------------------------------------------%                                                                           
load('oldresultsasa1.mat')                                                                           
subresults=oldresultsasa(oldresultsasa(:,2)<quantile(oldresultsasa(:,2),0.25),:);    
                                                                           %[Take the 25th percentile above the minimum found, and use it as the 
                                                                           % initial threshold, and take all the parmeters for which the objective 
                                                                           % function is below this threshold ]
ncol=size(subresults,1);
ciparameter=subresults(randperm(ncol,200),5:end);                          %[Draw a random sample from all selected parameters
                                                                           % size=200]

save('ciparameter','ciparameter');                                         %
clear ncol                                                                 %
                
                                                                           
subsamples(datairline)                                                     %[Draw 100 random samples from the orginial dataset, 
                                                                           % size=1/4 original size]

                                                                           
replication_confidencefcn();                                               %[Get the threshold for confidence region, one can redo step 
                                                                           % 12 for several times, until there is no big changes in c]

                                                                           
load savedvalues1;                                                                                                                                                  load savedvalues1;
temp(1,:)=[];                                                              % Drop first row and first column which are just counting.
temp(:,1)=[];                                                              %
                                                                           %
maxofrow=zeros(size(temp,1),1);                                            %
maxofcol=zeros(1,size(temp,2));                                            %
                                                                           %
colquant=zeros(1,size(temp,2));                                            %
rowquant=zeros(size(temp,1),1);                                            %
                                                                           %
for i=1:size(temp,2)                      
   
    sortbycol=sort(temp(:,i));                                             %
    colquant(i)=sortbycol(floor(size(temp,1)*0.95));                       % Choose the 95 percent threshold.
end                                                                        %
                                                                           %
bycolcutoffpnt=max(colquant);                                              % This is the desired outcome.
                                                                           %
                                                                                                                 
CR=oldresultsasa(oldresultsasa(:,2)<fval+bycolcutoffpnt/rowX,:);
                                                                           
%--------------------------------------------------------------------------%
% (13) CONFIDENCE INTEVAL -- PROJECT CR ONTO EACH COVARIATES               %
%--------------------------------------------------------------------------% 

LB=min(CR(2:end,5:end));
UB=max(CR(2:end,5:end));
