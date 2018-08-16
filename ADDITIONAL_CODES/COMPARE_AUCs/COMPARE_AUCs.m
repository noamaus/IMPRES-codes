load('tab2.mat')
load('SIZES.mat')
load('pred.mat')
T=tab2;

%%get sizes for weighted means
wights = [length(SIZES.P1)+length(SIZES.N1),length(SIZES.P2)+length(SIZES.N2),length(SIZES.P3)+length(SIZES.N3),length(SIZES.P4)+length(SIZES.N4),length(SIZES.P5)+length(SIZES.N5),length(SIZES.P6)+length(SIZES.N6),length(SIZES.P7)+length(SIZES.N7),length(SIZES.P10)+length(SIZES.N10),length(SIZES.P11)+length(SIZES.N11)]';

prep = [2,5,7,8];%%pre PD-1
onp = [4,9];%% on PD-1
pd = [2,4,5,7,8,9];%% all PD-1
ct = [1,3,6];%%CTLA-4
MM=[];
SD=[];

%%get mean and std of AUCs per each group
for i = 1:8
    v=T(i,:);
    MM(i,4) = v(ct)*wights(ct)/sum(wights(ct));
    SD(i,4) = std(v(ct));
    MM(i,2) = v(prep)*wights(prep)/sum(wights(prep));
    SD(i,2) = std(v(prep));
    MM(i,3) = v(onp)*wights(onp)/sum(wights(onp));
    SD(i,3) = std(v(onp));
    MM(i,1) = v*wights/sum(wights);
    SD(i,1) = std(v);
    
    %%for overlaid scatters 
    SP1(i,:) = v(ct);
    SP2(i,:) = v(onp);
    SP3(i,:) = v(prep);
    SP4(i,:) = v;
end

%%compare performances all treatments (P-values on top of bars)
pp1=[]; 
for i = 1:8
    v1=T(i,:);
    v2=T(8,:);
    [pp1(i),hh] = signrank(v1,v2,'tail','left');
end

% pp1(pp1==pp1(1))=0.004;
[sv,si] = sort(MM(:,1));
SD1 = SD(si,:);

MM1 = MM(si,:);
pp2 = pp1(si);

bar(1:8,MM1,0.8)
ylabel('IMPRES')
% text(1:8,MM1,pred(si), 'HorizontalAlignment','center','VerticalAlignment','bottom')

pvt={};
for i = 1:8
    pvt{i} = ['P-value = ',num2str(pp2(i))];
end
z=[0.7:1:7.7;0.9:1:7.9;1.1:1:8.1;1.3:1:8.3];
hold on
errorbar(z',MM1,SD1,'.')

text(1:8,sv,pvt,'VerticalAlignment','top')


SP11 = SP1(si,:);
SP22 = SP2(si,:);
SP33 = SP3(si,:);
SP44 = SP4(si,:);
xa1=reshape(repmat(z(4,:),3,1),24,1);
ya1=reshape(SP11',24,1);
scatter(xa1,ya1,20,'filled','MarkerFaceColor',[0 .9 .7]);

xa1=reshape(repmat(z(2,:),4,1),32,1);
ya1=reshape(SP33',32,1);
scatter(xa1,ya1,20,'filled','MarkerFaceColor',[0.7 .7 .1]);


xa1=reshape(repmat(z(3,:),2,1),16,1);
ya1=reshape(SP22',16,1);
scatter(xa1,ya1,20,'filled','MarkerFaceColor',[0 .1 .6]);



xa1=reshape(repmat(z(1,:),9,1),72,1);
ya1=reshape(SP44',72,1);
scatter(xa1,ya1,20,'filled','MarkerFaceColor',[0.8 .1 .1]);

