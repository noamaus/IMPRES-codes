function [IMPRES,signd3] = getIMPRESRAT(Data,CPall,F)

%%%GET IMPRES USING FEATURE RATIOS INSTEAD OF BINARY RELATIONS

AUC = [];

%%%Create table of comparisons
[a,b,c]=intersect(CPall,Data.genes);
for k = 1:length(Data.sample)
    cnt=1;
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(j);
            g2(cnt) = a(i);
            signd3(cnt,k) = Data.GE(c(i),k)/Data.GE(c(j),k);
            cnt = cnt+1;
        end
    end
    
    %%%
    for i = 1:(length(c)-1)
        for j = i+1:length(c)
            g1(cnt) = a(i);
            g2(cnt) = a(j);
            signd3(cnt,k) = Data.GE(c(j),k)/Data.GE(c(i),k);
            cnt = cnt+1;
        end
    end
    %%%
end
[aa,bb] = size(signd3);

IMPRES = (signd3(F,:));

