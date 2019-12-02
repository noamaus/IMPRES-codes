function [FEATS] = SelectFeaturs(RES,CFF)

T1 = RES.FEATS(:,RES.AUCTS>=0.6);

scorePos = sum(T1');


T2 = RES.FEATS(:,RES.AUCTS<=0.4);
scoreNeg = sum(T2');
scoreTot = scorePos-scoreNeg;

quant = quantile(scoreTot(CFF),0.95);
FEATS = CFF(scoreTot(CFF)>=quant);

end

