%%%PART1: perform feature selection 
%%%to define IMPRES score **NO MELANOMA DATA**


%[IMPRES, RAT] = getIMPRES(NEUB,CPall,1:756) - get initial relations tab in NB for all considered
%features

%[RES] = hilClibming(RATS,FOLDS.STSP,FOLDS.STSN,FOLDS.STP,FOLDS.STN,CFF,NEUB,CPall)
%perform the hill clibming feature section procedure via 500 repetitions of
%5-fold cross validation with folds saved in FOLDS
%***REQUIRES LONG RUNNING TIME***4-6 hours

%[FEATS] = SelectFeaturs(RES,CFF); - select the optimal features resulting
%from the feature selection procedure (top 0.05 scored features, binomial p<=0.05)

