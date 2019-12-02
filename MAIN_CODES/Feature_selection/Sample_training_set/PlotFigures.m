%%%%%%%%A. PLOT_FIGURES:
%%Plotting (1) the fraction of features recurring in the 500 feature sets, and 
%%(2) Boxplots showing the area under the ROC curve (AUC) for responder/non-responder 
%%classification obtained with these 200 features sets selected over the NB data, 
%%for each of the melanoma datasets



load('CFF.mat')
load('RR.mat') %%load 500 sampled feature sets

% %%PART3: evaluate melanoma performance 
for i = 1:length(RR.RES)
    i
    [FEATS] = SelectFeaturs(RR.RES(i),CFF);
    [AUCall(:,i)] = PredictAllMelanomaData(FEATS);
end



%%%%If more than one training is sampled:
load('CNTFT.mat')
if length(RR.RES)>1
    subplot(5,1,1:2)
    cntF = sum(CNTFT')/500;
    [sv1,si1] = sort(cntF,'descend');
    bar(cntF(si1(1:30)))
    ylabel('Fraction corrected feature sets')

    subplot(5,1,3:5)
    load('D.mat')
    [sv,si] = sort(median(AUCall'),'descend');
    boxplot(AUCall(si,:)',D(si),'LabelOrientation','inline')
    ylim([0,1.1])
    ylabel('AUC')

end