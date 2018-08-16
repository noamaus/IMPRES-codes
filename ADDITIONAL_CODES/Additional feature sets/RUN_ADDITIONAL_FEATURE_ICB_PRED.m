load('A_FEAT.mat')%%%other feature sets


%%% predict ICB response for each feature set
subplot(2,2,1)
AUCall1 = PredictICBresponse(A_FEAT.f1);
title('selected for Hugo and Van Allen et al')

subplot(2,2,2)
AUCall1 = PredictICBresponse(A_FEAT.f2);
title('selected for Hugo and Riaz et al')

subplot(2,2,3)
AUCall1 = PredictICBresponse(A_FEAT.f3);
title('selected using both PCA clusters')


subplot(2,2,4)
AUCall1 = PredictICBresponse(A_FEAT.f4);
title('reduced feature sets')
