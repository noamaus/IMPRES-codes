%%%%%%%% Sample NB controlled, unbiased training set (as of IMPRES
%%%%%%%% training set)


load('POOL.mat')
load('CFF.mat')
load('NEUB.mat')
load('CPall.mat')

s1 = sum(POOL.FIRSTS(:,POOL.AUCTS>=0.6)'); %%Positive draws
s2 = sum(POOL.FIRSTS(:,POOL.AUCTS<=0.4)'); %%Negative draws
ss = sum(POOL.FIRSTS');

fp = find(s1>0); %%initial feature of folds
fn = find(s2>0);
f11 = find(ss);

sfp=s1(f11);
sfn=s2(f11);
sfA=ss(f11);


% % % FOR RANDOM PARTITION

rounds = 1; %% Change if want to derive more than one training set

% rounds = 200; 

CNTFT=zeros(756,rounds);


l2 = length(POOL.AUCTS);
for i = 1:rounds
    
    DRS = ones(size(f11));
    NPSd=zeros(size(f11));
    NNSd=zeros(size(f11));
    
    RN1 = rng; %%%To save random generator
    Pdraws = 100+randi(200);
    rng(RN1);
    Ndraws = randi(500-1-Pdraws);
    
    
    %%%GET RANDOM DISTRIBUTION FOR SAMPLING:
    RTI=50; %%%CAN CHANGE - features often initiating (may sample more by random)
    
    
    RN1 = rng; %% TO Save randon generator
    
    
    rng(RN1);
    TH1 = randi(floor(min(quantile(sfA,0.75),RTI)));
    RLI = find(sfp>TH1&sfn>TH1); %%%CHANGE
    rng(RN1);
    RA1 = randi(length(RLI));
    rng(RN1);
    rp1=randperm(length(RLI),RA1);
    RLI2 = RLI(rp1);
    
    
    
    
    rng(RN1);
    NPSd(RLI2) = randi(TH1,size(RLI2));
    rng(RN1);
    NNSd(RLI2) = randi(TH1,size(RLI2));
    
    
    while sum(NPSd)<Pdraws
        z = ceil(min(((Pdraws-sum(NPSd))/2),(nnz(NPSd<sfp)/2)));
        fk=find(NPSd<sfp);
        rng(RN1);
        rp=randperm(length(fk),z);
        NPSd(fk(rp)) = NPSd(fk(rp))+1;
    end

    NNSd(NPSd==0) = 1;
    while sum(NNSd)<Ndraws
        z = ceil(min(((Ndraws-sum(NNSd))/2),(nnz(NNSd<sfn)/2)));
        fk=find(NNSd<sfn);
        rng(RN1);
        rp=randperm(length(fk),z);
        NNSd(fk(rp)) = NNSd(fk(rp))+1;
    end

    
    
    
    i

    
    IPP = [];INN = [];
    for j = 1:length(f11)

        IDP = NPSd(j);%%number of positive draws with j as initial feature
        IDN = NNSd(j);%%number of negative draws with j as initial feature
        

        ip = (POOL.FIRSTS(f11(j),:)==1);%%%j initial feature
    
        ip1 = find((POOL.AUCTS>=0.6)&ip);%%%Positive indices with j initial feature
        in1 = find((POOL.AUCTS<=0.4)&ip);%%%Negative indices with j initial feature
        ia1 = find((POOL.AUCTS>0.4)&(POOL.AUCTS<0.6));%%%Neutral indices 
        
        
        %%Do the sampling of positive
        if IDP>0&nnz(ip1)>0
            rng(RN1);
            ri1 = randi(length(ip1),IDP,1)';
            IPP = [IPP,(ip1(ri1))];
        end
        
         %%Do the sampling of negative
        if IDN>0&nnz(in1)>0
            rng(RN1);
            ri1 = randi(length(in1),IDN,1)';
            INN = [INN,(in1(ri1))];    
        end
        
    
    end
    
    
    
    cnts = nnz([IPP,INN]);
    rng(RN1);
    rn = randi(length(ia1),500-cnts,1);
    IA = rn';%%% Set neutral draws - optional
   
      
    %%% build the struct of sampled training set
    [REST] = CUTRES(POOL,[IPP,INN,IA]); %%%Unbiased sampled training set
     
     
    [IMPRES,RATS] = getIMPRES(NEUB,CPall,1:756);
    
    %%%run hill clibming to get selected features from sampled training
    %%%set
    [REST] = hilClibming(RATS,REST.STSP,REST.STSN,REST.STP,REST.STN,CFF,NEUB,CPall);
%     %%TO VALIDATE
    
    
    RR.RES(i) = REST; %%%SAVE THE RESULT (if more than one iteration)
    
    [FEATST] = SelectFeaturs(REST,CFF); 
    
	CNTFT(FEATST,i) = 1; %%%SAVE THE SELECTED FEATURES (if more than one iteration)
    
    
end






    
