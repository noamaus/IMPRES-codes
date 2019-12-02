function [AUCall] = PredictAllMelanomaData(FEATS)


load('CPall.mat')
load('CFF.mat')




%% load all melanoma Datasets
toPlot=0;
load('SKCM.mat')
load('VANALLEN.mat')
load('WARGO.mat')
load('HUGO.mat')
load('FELIP.mat')
load('TIMCHAN.mat')
load('MGH.mat')




%%%1. Van Allen data aCTLA-4
hold on 
P1 = find(strcmp(VANALLEN.response_CTLA4,'response'));
N1 = find(strcmp(VANALLEN.response_CTLA4,'nonresponse')|strcmp(VANALLEN.response_CTLA4,'long-survival'));
[AUC1, score1, label1, jj] = classifyImmuneCOMP(VANALLEN,CPall,FEATS,P1,N1,toPlot);


%%%2. Hugo data aPD-1
hold on 
P2 = find(strcmp(HUGO.response_PD1,'Complete Response')|strcmp(HUGO.response_PD1,'Partial Response'));
N2 = find(strcmp(HUGO.response_PD1,'Progressive Disease'));
[AUC2, score2, label2, jj] = classifyImmuneCOMP(HUGO,CPall,FEATS,P2,N2,toPlot);

%%%TCGA SKCM data aCTLA-4
hold on 
P3 = find(strcmp(SKCM.response,'Complete Response')|strcmp(SKCM.response,'Partial Response'));
N3 = find(strcmp(SKCM.response,'Clinical Progressive Disease')|strcmp(SKCM.response,'Stable Disease'));
[AUC3, score3, label3, jj] = classifyImmuneCOMP(SKCM,CPall,FEATS,P3,N3,toPlot);


%%%Wargo data on aPD-1 (post aCTLA4)
hold on 
P4 = find(strcmp(WARGO.antiPD1,'R')&(~cellfun(@isempty, strfind(WARGO.time,'On aPD1'))));
N4 = find(strcmp(WARGO.antiPD1,'NR')&(~cellfun(@isempty, strfind(WARGO.time,'On aPD1'))));
[AUC4, score4, label4, jj] = classifyImmuneCOMP(WARGO,CPall,FEATS,P4,N4,toPlot);


%%%Wargo data pre aPD-1 (post aCTLA4)
hold on 
P5 = find(strcmp(WARGO.antiPD1,'R')&(~cellfun(@isempty, strfind(WARGO.time,'Pre aPD1'))));
N5 = find(strcmp(WARGO.antiPD1,'NR')&(~cellfun(@isempty, strfind(WARGO.time,'Pre aPD1'))));
[AUC5, score5, label5, jj] = classifyImmuneCOMP(WARGO,CPall,FEATS,P5,N5,toPlot);
% score5 = (15/14)*score5;


%%%Wargo data pre aCTLA4
hold on 
P6 = find(strcmp(WARGO.antiCTLA4,'R')&strcmp(WARGO.time,'Pre-aCTLA4'));
N6 = find(strcmp(WARGO.antiCTLA4,'NR')&strcmp(WARGO.time,'Pre-aCTLA4'));
[AUC6, score6, label6, jj] = classifyImmuneCOMP(WARGO,CPall,FEATS,P6,N6,toPlot);


%%%%FELIP pre aPD-1
hold on 
P7 = find(strcmp(FELIP.response,'CR')|strcmp(FELIP.response,'PR'));
N7 = find(strcmp(FELIP.response,'PD'));
[AUC7, score7, label7, jj] = classifyImmuneCOMP(FELIP,CPall,FEATS,P7,N7,toPlot);

%%%%TIM CHAN pre a-PD-1
hold on 
P10 = find((strcmp(TIMCHAN.response,'PR')|strcmp(TIMCHAN.response,'CR'))&(~cellfun(@isempty, strfind(TIMCHAN.sample', 'Pre'))));
N10 = find((strcmp(TIMCHAN.response,'PD')|strcmp(TIMCHAN.response,'SD'))&(~cellfun(@isempty, strfind(TIMCHAN.sample', 'Pre'))));
[AUC10, score10, label10, jj] = classifyImmuneCOMP(TIMCHAN,CPall,FEATS,P10,N10,toPlot);

%%%%TIM CHAN on a-PD-1
hold on 
P11 = find((strcmp(TIMCHAN.response,'PR')|strcmp(TIMCHAN.response,'CR'))&(~cellfun(@isempty, strfind(TIMCHAN.sample', 'On'))));
N11 = find((strcmp(TIMCHAN.response,'PD')|strcmp(TIMCHAN.response,'SD'))&(~cellfun(@isempty, strfind(TIMCHAN.sample', 'On'))));
[AUC11, score11, label11, jj] = classifyImmuneCOMP(TIMCHAN,CPall,FEATS,P11,N11,toPlot);


hold on 
P12 = find(strcmp(MGH.response_CTLA4,'responder')|strcmp(MGH.response_CTLA4,'PR')|strcmp(MGH.response_CTLA4,'SD'));
N12 = find(strcmp(MGH.response_CTLA4,'mixed / PD')|strcmp(MGH.response_CTLA4,'PD')|strcmp(MGH.response_CTLA4,'SD / no response')|strcmp(MGH.response_CTLA4,'non-responder')|strcmp(MGH.response_CTLA4,'none'));
[AUC12, score12, label12, jj] = classifyImmuneCOMP(MGH,CPall,FEATS,P12,N12,toPlot);

hold on 
P13 = find(strcmp(MGH.response_PD1,'responder')|strcmp(MGH.response_PD1,'PR')|strcmp(MGH.response_PD1,'SD'));
N13 = find(strcmp(MGH.response_PD1,'mixed / PD')|strcmp(MGH.response_PD1,'PD')|strcmp(MGH.response_PD1,'SD / no response')|strcmp(MGH.response_PD1,'non-responder')|strcmp(MGH.response_PD1,'none'));
[AUC13, score13, label13, jj] = classifyImmuneCOMP(MGH,CPall,FEATS,P13,N13,toPlot);




AUCall = [AUC1;AUC2;AUC3;AUC4;AUC5;AUC6;AUC7;AUC10;AUC11;AUC12;AUC13];

