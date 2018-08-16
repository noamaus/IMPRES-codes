function [AUC, score, label, XX] = classifyImmuneCOMP(Data,CPall,F,Pos,Neg,plotROC)

%%%Classify based on comparisons of the F pairs from CPall generated table        
%%%Pos are responders to immunotherapy, Neg are non-responders
AUC = [];

%%%Create table of comparisons
[a,b,c]=intersect(CPall,Data.genes);
for k = 1:length(Data.sample)
    cnt=1;
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(j);
            g2(cnt) = a(i);
            signd3(cnt,k) = Data.GE(c(i),k)<Data.GE(c(j),k);
            diff3(cnt,k) = Data.GE(c(i),k)-Data.GE(c(j),k);
            cnt = cnt+1;
        end
    end
    
    %%%
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(i);
            g2(cnt) = a(j);
            signd3(cnt,k) = Data.GE(c(i),k)>Data.GE(c(j),k);
            diff3(cnt,k) = Data.GE(c(j),k)-Data.GE(c(i),k);
            cnt = cnt+1;
        end
    end
    %%%
end

X0 = signd3(:,Pos);
X1 = signd3(:,Neg);

Y0 = diff3(:,Pos);
Y1 = diff3(:,Neg);

xx0 = X0(F,:);
xx1 = X1(F,:);

yy0 = Y0(F,:);
yy1 = Y1(F,:);

XX = [xx0,xx1];
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


