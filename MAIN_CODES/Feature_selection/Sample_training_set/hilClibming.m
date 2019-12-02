function [RES] = hilClibming(RATS,STSP,STSN,STP,STN,CFF,NEUB,CPall)

%%%Runs Hill climbming feature selection with five fold cross validation
%%%for each fold (up to 15 steps each fold)

%Input:
%%%RATS: tables of all considered checkpoint genes ratios caculated before
%%%STSP,STSN,STP,STN: are the 500 saved folds (TSTP/N are training sets and
%%%STN/P are the saved test sets)
%%%CFF: are condifered pairs (with one checkpoint genes related to
%%%CTLA-4/PD-1 in the pair)
%%%NEUB: the neuroblastoma data (normalized and filtered)
%%%CPall: checkpint genes

%Output:
%%%RES: a structure containing saved training/test sets, AUCTS: the
%%%resutling AUCs from 500 folds on the training sets, FEATS: the selected
%%%features for each of the 500 folds

reps = 500;
FEATS = zeros(756,reps);

for k = 1:reps
    k

   
    TSN = STSN(:,k);
    TSP = STSP(:,k);
    TN = STN(:,k);
    TP = STP(:,k);
    
    
    
    curAUC = 0.1;
    mAUC = 0;
    round = 1;
    csel=[];

    while round<=15

        for i = 1:length(CFF)
            c2 = unique([csel; CFF(i)]);

            [AUCTR(i)] = classifyImmuneCOMP3(NEUB,CPall,c2,TSP,TSN,0,RATS);
        end
        ni = find(AUCTR==max(AUCTR));
        curAUC = max(AUCTR);
        if mAUC<curAUC
            csel = [csel;CFF(ni)];
            round = round+1;
            mAUC = curAUC;
        else
            break
        end
    end
    [AUCTS(k)] = classifyImmuneCOMP3(NEUB,CPall,c2,TP,TN,0,RATS);
    FEATS(csel,k) = 1;
    STSN(:,k) = TSN;
    STSP(:,k) = TSP;
    STN(:,k) = TN;
    STP(:,k) = TP;
end

RES.STSP=STSP;
RES.STSN=STSN;
RES.STP=STP;
RES.STN=STN;
RES.AUCTS = AUCTS;
RES.FEATS=FEATS;


    

