load('CPall.mat') %%%checkpint genes
load('NEUB.mat') %%Neublastoma data 

%%% Obtain immune relation tab for NB data (all features)
[IMPRES,RATS] = getIMPRES(NEUB,CPall,1:756);

load('CFF.mat') %%%considered relations (containing CTLA-4/PD-1 related genes as at least one gene in the pair)
load('FOLDS.mat') %%%the folds structure to perform hill clibming feature selection on NB data 


%%%***Warning: moderate-long running time (4-6 hours) (result is saved as a structure RES.mat)***
% [RES] = hilClibming(RATS,FOLDS.STSP,FOLDS.STSN,FOLDS.STP,FOLDS.STN,CFF,NEUB,CPall)


load('RES.mat') %%%output from 'hilClibming.m' 
% % %%%select best features from previous procedure
[FEATS] = SelectFeaturs(RES,CFF);


%%% ROC(and AUC) for the selected  features with the full NB data
P8 = find(NEUB.high_risk==0&NEUB.age<=1.5&NEUB.progression==0);
N8 = find(NEUB.high_risk==1&NEUB.age<=1.5&NEUB.progression==1);
[AUCNB, scoreNB, labelNB, jj] = classifyImmuneCOMP(NEUB,CPall,FEATS,P8,N8,1);
