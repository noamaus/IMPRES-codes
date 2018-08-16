function [AUC, score, label] = classifyImmuneCOMP3(Data,CPall,F,Pos,Neg,plotROC,RatioTab)
%%%fast classification (when the relation tab is previously calculated)
%%%Classify based on comparisons of the F pairs (from the considered pairs in CPall)
%generated table of pairs   

%Input:
%%%Data: a structure with all the data
%%%CPall: are checkpoint genes
%%%F are the features to use for classification
%%%Pos: are cancer regression, Neg: are progressive disease
%%%plotROC: indicator if to plot ROC curve (1<-->plot)
%%%RatioTab: table with all possible ratios of considered checkpint genes in data

AUC = [];

[a,b,c]=intersect(CPall,Data.genes); 


X0 = RatioTab(F,Pos);
X1 = RatioTab(F,Neg);

XX = [X0,X1];

if length(F)>1
    score = sum(XX/length(F));
else
    score = double(XX);
    
end
  
label=[ones(1,length(Pos)),zeros(1,length(Neg))];

[XX1,YY1,T,AUC] = perfcurve(label',score',1);


if plotROC
    plot(XX1,YY1);
    xlabel('False positive rate')
    ylabel('True positive rate')
end


