%%%Get features from part 1
load('FEATS.mat')
load('CPall.mat')

%%%Load survival structures
load('TCH_surv.mat')
load('VAL_surv.mat')
load('HU_surv.mat')
load('FELIP_surv.mat')

%%%TIM CHAN data KM analysis 
subplot(2,4,1:2)
[IMPRES] = getIMPRES(TCH_surv,CPall,FEATS);
[p1, dAUC1] = SurvivalAnalysisForMELANOMA(TCH_surv,IMPRES)
title('Riaz et al')

%%%VAN ALLEN data KM analysis 
subplot(2,4,3:4)
[IMPRES] = getIMPRES(VAL_surv,CPall,FEATS);
[p2, dAUC2] = SurvivalAnalysisForMELANOMA(VAL_surv,IMPRES)
title('Van Allen et al')


%%%PFS comparison PRAT(FELIP) data (high/low IMPRES)
subplot(2,4,5)
[IMPRES] = getIMPRES(FELIP_surv,CPall,FEATS);
[p3,h3] = ranksum(FELIP_surv.PFS(IMPRES>=median(IMPRES)),FELIP_surv.PFS(IMPRES<median(IMPRES)),'tail','right')
boxplot(FELIP_surv.PFS(IMPRES>=median(IMPRES)))
title('Prat et al')
ylim([0,33])
xlabel('HIGH IMPRES')
ylabel('PFS (months)')
subplot(2,4,6)
boxplot(FELIP_surv.PFS(IMPRES<median(IMPRES)))
ylim([0,33])
title('Prat et al')
xlabel('LOW IMPRES')
ylabel('PFS (months)')

%%%PFS comparison VAN ALLEN data (high/low IMPRES)
subplot(2,4,7)
[IMPRES] = getIMPRES(VAL_surv,CPall,FEATS);
[p4,h4] = ranksum(VAL_surv.PFS(find(VAL_surv.IMPRES>median(VAL_surv.IMPRES))),VAL_surv.PFS(find(VAL_surv.IMPRES<=median(VAL_surv.IMPRES))),'tail','right')
ylim([0,50])
boxplot(VAL_surv.PFS(find(VAL_surv.IMPRES>median(VAL_surv.IMPRES))))
title('Van Allen et al')
xlabel('HIGH IMPRES')
ylabel('PFS (months)')
subplot(2,4,8)
boxplot(VAL_surv.PFS(find(VAL_surv.IMPRES<=median(VAL_surv.IMPRES))))
ylim([0,50])
title('Van Allen et al')
xlabel('LOW IMPRES')
ylabel('PFS (months)')



