function[p, dAUC] = SurvivalAnalysisForMELANOMA(Data,vals)
%%%Kaplan Meier survival analysis

%Input:
%%%Data - Data structure with fields: 
%%%survival:survival time (days), and death: 1 if censored, 0 if dead
%%%vals: values(score) to use to seperate survival groups

%Input:
%%%p: logrank p-value and dAUC: delta AUC resulting from KM survival
%%%analysis

[A,idx] = sort(vals);
low = quantile(vals,0.5);
high = quantile(vals,0.5);

if nnz(A<low)>nnz(A>high)
    a = find(A<low);
    b = find(A>=high);
else
    a = find(A<=low);
    b = find(A>high);
end

aa = size(a);
bb = size(b);
Data.survival = Data.survival(idx);
Data.death = Data.death(idx);


[p]=logrank([Data.survival(1:a(end)), Data.death(1:a(end))],[Data.survival(b(1):end), Data.death(b(1):end)],0.5);

m = mean(Data.survival(1:a(end)))>mean(Data.survival(b(1):end));



[t1,T1]=kmplot([Data.survival(1:a(end)) Data.death(1:a(end))]);
S1=stairs(t1,T1,'b'); 
xmax=max(t1)+1;
hold on

[t2,T2]=kmplot([Data.survival(b(1):end) Data.death(b(1):end)]);
S1=stairs(t2,T2,'m');
xmax=max(t2)+1;


totA = max(max(t1),max(t2));
AUC1 = MyAUC(t1,T1)/totA;
AUC2 = MyAUC(t2,T2)/totA;
dAUC = (AUC1-AUC2);

xlabel('Followup (Days)');
ylabel('Overall survival');
legend('Low IMPRES','High IMPRES');



