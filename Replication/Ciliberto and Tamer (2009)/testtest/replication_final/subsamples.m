% SUBSAMPLES randomly take parts of the original dataset, and meanwhile
% collect the corresponding empirical probabilities and unobervables for
% each observations. The size of subsample is 1/4 of the original one.
% Here, for time's sake, I only draw 25 subsamples.

function subsamples(datairline)

global rowX prob epsi

rng(4)
rndcol=randperm(rowX,floor(rowX/4));
subsample=datairline(rndcol,:);
subprob=prob(rndcol,:);
subepsi=epsi(rndcol,:);
clear rndcol;

for i=2:25
    
rng(3+i);
rndcol=randperm(rowX,floor(rowX/4));
subsample=[subsample;datairline(rndcol,:)];
subprob=[subprob; prob(rndcol,:)];
subepsi=[subepsi; epsi(rndcol,:)];
clear rndcol;
    
end

save('subsample','subsample');
save('subprob','subprob');
save('subepsi','subepsi');

