%%%Survival analysis: Perform survival analysis using IMPRES score
%%%(presented in Figure 2)

%%%Loading survival structure (structures containing expression, survival
%%%time and death for samples with survival imformation)


%%%[IMPRES] = getIMPRES(Data_surv,CPall,FEATS): Calculate IMPRES for
%%%samples with survival info using FEATS - the feature selected in PART1
%%%and CPall: immune checkpoint genes

%%%[p, dAUC] = SurvivalAnalysisForMELANOMA(Data,IMPRES): PErform KM survival
%%%analysis using IMPRES score for groups defenition
