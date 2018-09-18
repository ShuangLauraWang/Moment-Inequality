%Edited by Shuang Wang for EC711 Final Project
%March 18, 2017


%--------------------------------------------------------------------------%
% CONFIDENCEFCN computes the treshold needed to construct the confidence   %
% intervals                                                                %
%                                                                          %
%   function confidencefcn()                                               %
%                                                                          %
%                                                                          %
% NOTES                                                                    %
%       (i) This is written as a function because we later on compile the  %
%       files in an unique executable. However, this can be run            %
%       interactively within matlab. In that case the first line, function %
%       mainhete(), can be commented out.                                  %
%       (ii) It is written in a way that makes it easy to run many of      %
%       these in parallel without having to use paralle programming.       %
%       (iii) There are no outputs or inputs. The file uploads data and    %
%       save the results in a matlab file.                                 %
%       (iv) We only comment the lines that are new relative to the file   %
%       MAINHETE.m                                                         %
%                                                                          %
%                                                                          %
% Written by Federico Ciliberto and Elie Tamer, March 2003.                %
% When using this code or parts of it, please cite the following article:  %
% Federico Ciliberto and Elie Tamer, "Market Structure and Multiple        %
%          Equilibria in the Airline Industry," Econometrica, Vol. 77,     %
%          No. 6 (November, 2009), 1791-1828.                              %
%                                                                          %
%--------------------------------------------------------------------------%
                                                                           %
function replication_confidencefcn()                                                   %
                                                                           %
%--------------------------------------------------------------------------%
% (1) READS IN THE INITIAL VALUES OF THE PARAMETERS AND AN INDICATOR THAT  %
% CAN BE ASSOCIATED TO THEM                                                %
%--------------------------------------------------------------------------%
                                                                           %
start_values=textread('start_values');                                     %
                                                                           %
sim_num=num2str(start_values(1,1));                                        %
loopstart=start_values(1,1);                                               % This line is used to have parallel computers working independently.
                                                                           %
param0 = [start_values(2,1) start_values(2,10:17)]';  
                                                                           %
load ('subsample.mat')
load('subprob.mat')
load('subepsi.mat')                                                        %
%--------------------------------------------------------------------------%
% (2) PREPARE THE MATRIX WHERE THE RESULTS ARE SAVED                       %
%--------------------------------------------------------------------------%
                                                                           %
temp=[0,zeros(1,200)];                                                     % This is the matrix that we need to use to get the thresholds.
                                                                           % Here we are looking at 100 subsamples and 200 parameters.
                                                                           % Dimension: # subsamples (e.g. 100) X # of (sample of) parameters satisfying the 
                                                                           %                                       conjectured treshold condition (e.g. 200)
save(['savedvalues',sim_num],'temp');                                      % This line saves the matrix above.
clear temp                                                                 %
                                                                           %
for m=(1+25*(loopstart-1)):(25*loopstart)                                  % 25 says that 4 computers were running in parallel,
                                                                           % each dealing with 25 subsamples.
                                                                           %
%--------------------------------------------------------------------------%
% (3) LOADS DATA                                                           %
%--------------------------------------------------------------------------%
    mstr=num2str(m);                                                       % This line records which subsample we are looking at.
    global subsize
    
    subsize=size(subsample,1)/25;
    datairline=subsample(1+(m-1)*subsize:subsize+(m-1)*subsize,:);       %
                                                                           %    
%--------------------------------------------------------------------------%
% (4) DEFINE VARIABLE COUNTING THE NUMBER OF FIRMS, THE MATRIX OF POSSIBLE %
% MARKET STRUCTURES                                                        %
%--------------------------------------------------------------------------%
                                                                           %
    global index total k                                                   %
    k=2;                                                                   %
    index=makeindex(k);                                                    %
    total=size(index,1);                                                   %
                                                                           %
%--------------------------------------------------------------------------%
% (5) SET THE NUMBER OF SIMULATIONS                                        %
%--------------------------------------------------------------------------%
                                                                           %
    global r                                                               %
    r=100;                                                                 %
                                                                           %
%--------------------------------------------------------------------------%
% (6) DEFINE VARIABLES ASSOCIATED WITH DATA                                %
%--------------------------------------------------------------------------%
                                                                           %
global X Y rowX otherentry                                                 %

Xheter=double(datairline(:,{'marketpresenceAA','marketpresenceDL'...
    'mindistancefromhubAA', 'mindistancefromhubDL'}));                     % Exogenous regressors.

otherentry=double(datairline(:,{'airlineUA','airlineAL','airlineLCC',...
    'airlineWN'}));                                                        % The entry decisions of UA, MA, LCC & WN 
                                                                           % are taken as given.

X=[Xheter otherentry];                                                     % Dimension: # markets X # regressors
                                                                           %
colX=size(X,2);                                                            % Number of exogenous regressors.
                                                                           %
Y=double(datairline(:,{'airlineAA','airlineDL'})); 

rowX=size(X,1); 
                                                   %
                                                                           %
%--------------------------------------------------------------------------%
% (7) FIRST STAGE EMPIRICAL PROBABILITIES FOR EACH SUBSAMPLE               %
%--------------------------------------------------------------------------%
                                                                           %
    global prob                                                            %
    prob=subprob(1+(m-1)*subsize:subsize+(m-1)*subsize,:);         %
                                                                           %
%--------------------------------------------------------------------------%
% (8) LOAD THE SIMULATED UNOBSERVABLES                                     %
%--------------------------------------------------------------------------%
    global epsi                                                            %
    epsi=subepsi(1+(m-1)*subsize:subsize+(m-1)*subsize,:);                                            %
                                                                           %
%--------------------------------------------------------------------------%
% (9) CONSTRUCTS INDEXES THAT PERMIT TO AVOID HAVING TO LOOP OVER MARKETS  %
%--------------------------------------------------------------------------%
                                                                           %
    global indexheter                                                      %
    tempheter=reshape(1:k*size(X,1),size(X,1),k);                          %
    indexheter=zeros(size(X,1)*total,k);                                   %
    for j=1:size(X,1)                                                      %
        indexheter((j-1)*total+1:(j-1)*total+total,:)=repmat(tempheter(j,:),total,1);
    end                                                                    %
    clear tempheter                                                        %
                                                                           %
    global repindex                                                        %
    repindex=repmat(index,size(X,1),1);                                    %
                                                                           %
    global diffdummies                                                     %
    diffdummies=repmat(reshape(1:total*k,total,k),size(X,1),1);            %
                                                                           %
%--------------------------------------------------------------------------%
% (10) LOAD THE SET OF PARAMETERS THAT WERE CHOSEN USING THE CONJECTURED   %
% THRESHOLD, AND COMPUTE THE DISTANCE FUNCTION FOR EACH ONE OF THESE       %
% PARAMETERS IN THE SUBSAMPLE UNDER CONSIDERATION                          %
%--------------------------------------------------------------------------%
                                                                           %
    load ciparameter;                                                      %
    cutoffcoeffs=200;                                                      %
                                                                           %
    compare=zeros(cutoffcoeffs,1);                                         %
    for i=1:200                                                            %
        [compare(i)]=replication_simuconf(ciparameter(i,:)');              % This has a different name but is basically the same function as simuhete.m
                                                                           % The only difference is that we do not save the parameters while minimizing
    end                                                                    %
                                                                           %
%--------------------------------------------------------------------------%
% (11) FIND THE ARGMIN OF THE DISTANCE FUNCTION FOR THIS SUBSAMPLE         %
%--------------------------------------------------------------------------%
                                                                           %
    lb = ones(1,9)*(-30);                                                  %
    ub = ones(1,9)*30;                                                     %
                                                                           %
    [param0,fval] = simulannealbnd(@replication_simuconf, param0 ,lb, ub); % Here a time limit might be necessary
                                                                           %    
%--------------------------------------------------------------------------%
% (12) FOR EACH OF THE PARAMETERS LOADED IN SECTION (9), COMPUTE THE       %
% DIFFERENCE BETWEEN THE VALUE THAT THE DISTANCE FUNCTION TAKE AT THOSE    %
% PARAMETERS AND THE MINIMIZED VALUE FOUND ABOVE                           %
%--------------------------------------------------------------------------%
                                                                           %
    fval=repmat(fval,cutoffcoeffs,1);                                      %
                                                                           %
    load(['savedvalues',sim_num])                                          %
                                                                           %
    temp=[temp;[m,(((compare-fval)/10000)*size(datairline,1))']];          %
                                                                           %
    save(['savedvalues',sim_num],'temp');                                  %
                                                                           %    
    clear temp                                                             %
    clear AirRndErr epsiairport epsimarket prob OPTIONS                    %
    clear fval r X i rowX coeffconfidsimple index compare                  %
    clear iteration start_values cutoffcoeffs k total                      %
    clear datairline epsi param                                            %
end                                                                        %
                                                                           %
end                                                                        %
                                                                           %
                                                                           %
                                                                           %                                                                           
                                                                           %                                                                           
                                                                           %                                                                           
                                                                           %                                                                           
%--------------------------------------------------------------------------%
% (13) HERE WE COMPUTE THE THRESHOLD TO BE APPROPRIATELY ADDED TO THE      %
% MINIMUM VALUE OF THE DISTANCE FUNCTION.                                  %
% WE NORMALLY RUN THIS LAST PIECE SEPARATELY. IT IS INCLUDED HERE FOR      %
% SIMPLICITY OF EXPOSITION                                                 %
%--------------------------------------------------------------------------%
%                                                                            %
% load savedvalues1;
% temp(1,:)=[];                                                              % Drop first row and first column which are just counting.
% temp(:,1)=[];                                                              %
%                                                                            %
% maxofrow=zeros(size(temp,1),1);                                            %
% maxofcol=zeros(1,size(temp,2));                                            %
%                                                                            %
% colquant=zeros(1,size(temp,2));                                            %
% rowquant=zeros(size(temp,1),1);                                            %
%                                                                            %
% for i=1:size(temp,2)                      
%     %
%     sortbycol=sort(temp(:,i));                                             %
%     colquant(i)=sortbycol(floor(size(temp,1)*0.95));                       % Choose the 95 percent threshold.
% end                                                                        %
%                                                                            %
% bycolcutoffpnt=max(colquant);                                              % This is the desired outcome.
%                                                                            %
%                                                                            