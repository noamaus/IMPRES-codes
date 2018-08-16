load('PAN_TCGA.mat')
figure


%%%Get mean IMPRES,mutation count
uc = unique(PAN_TCGA.canc);
MM=[];SDM=[];MI=[];SDI=[];
for i = 1:length(uc)
    i
    ids = find(strcmp(PAN_TCGA.canc,uc(i)));
    cnts(i) = nnz(ids);
    MM(i) = mean(PAN_TCGA.MUT(ids));
    SDM(i) = std(PAN_TCGA.MUT(ids));
    MI(i) = mean(PAN_TCGA.IMPRES(ids));
    SDI(i) = std(PAN_TCGA.IMPRES(ids));
end


%%%sort by IMPRES
[MI1,si] = sort(MI);
SDI1 = SDI(si);
MM1 = MM(si);
SDM1 = SDM(si);
cnts2 = cnts(si);

for i = 1:length(uc)
    ntx{i} = ['n=',num2str(cnts2(i))];
end

%%%Plot bar of IMPRES
subplot(8,1,1)
hold on
bar(1:23,MI1,0.8)
ylim([2,16])
errorbar(1:23,MI1,SDI1,'.')
% text(1:23,MI1,ntx,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom') 

[dotxI,dotyI] = GENdotPlot(PAN_TCGA.IMPRES,PAN_TCGA.canc,uc(si));
scatter(dotxI,dotyI,10,'filled')
ylabel('IMPRES')
% xlabel(uc(si))



%%%Plot bar of mutation count
subplot(8,1,2)
hold on
bar(1:23,MM1,0.8)
ylim([0,2000])

errorbar(1:23,MM1,SDM1,'.')
% text(1:23,MI1,ntx,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom') 

[dotxM,dotyM] = GENdotPlot(PAN_TCGA.MUT,PAN_TCGA.canc,uc(si));
hold on
scatter(dotxM,dotyM,10,'filled')
ylabel('Mutation count')
% xlabel(uc(si))

%%%Plot corr plor
subplot(8,1,3:8)
scatter(MI1,MM1,cnts(si),'filled','MarkerEdgeColor',[0 .1 .1],'MarkerFaceColor',[0 .9 .7])
[rho1,pp1]=corr(MI1',MM1','type','Spearman')
xlabel('IMPRES')
ylabel('Mutation count')
dx = 0.1; dy = 0.1;


for i = 1:length(uc)
    txts{i} = [uc(si(i));'';ntx(i)];
end
text(MI1-dx, MM1-dy, txts);






