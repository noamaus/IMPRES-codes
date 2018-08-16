function [AUC, score, label, jj] = classifyImmuneCOMP(Data,CPall,F,Pos,Neg,plotROC)

%%%classification with computing of relation tab
%%%Classify based on comparisons of the F pairs (from the considered pairs in CPall)
%%%generated table of pairs   

%Input:
%%%Data: a structure with all the data
%%%CPall: are checkpoint genes
%%%F are the features to use for classification
%%%Pos: are cancer regression, Neg: are progressive disease
%%%plotROC: indicator if to plot ROC curve (1<-->plot)
%%%RatioTab: table with all possible ratios of considered checkpint genes in data

%Input:
%%%AUC,score,label: resulting AUC, score and labels from classification
%%%jj: All pairs of checkpint genes


AUC = [];

%%%Create table of comparisons
[a,b,c]=intersect(CPall,Data.genes);
for k = 1:length(Data.sample)
    cnt=1;
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(j);
            g2(cnt) = a(i);
            RatioTab(cnt,k) = Data.GE(c(i),k)<Data.GE(c(j),k);
            cnt = cnt+1;
        end
    end
    
    %%%
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(i);
            g2(cnt) = a(j);
            RatioTab(cnt,k) = Data.GE(c(i),k)>Data.GE(c(j),k);
            cnt = cnt+1;
        end
    end
    %%%
end

X0 = RatioTab(F,Pos);
X1 = RatioTab(F,Neg);

XX = [X0,X1];

if length(F)>1
    score = sum(XX);
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


jj = [g1;g2]';
