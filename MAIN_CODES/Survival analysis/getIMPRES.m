function [IMPRES,RATS] = getIMPRES(Data,CPall,F)
%%%Get IMPRES and/or realtion tab (of relevant checkpoint genes) with no classification     

%Input:
%%%Data: a structure with all the data
%%%CPall: are checkpoint genes
%%%F are the features to used to construct IMPRES (or all features)

%Output:
%%%IMPRES: is the IMPRES score calculated, RATS: is the table of all ratios
%%%between pairs of checpoint genes (28 choos 2=756 pairs)
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

RATS = RatioTab;
IMPRES = sum(RatioTab(F,:))./length(F);
