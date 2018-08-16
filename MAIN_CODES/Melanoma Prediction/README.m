%%%PART2: running on melanoma data 
%%%using the IMPRES score defined in PART1 (for NB) as is
%%%to classify response/no response in melanoma data from ICB treated
%%%patients

%%%Loading the 6 sturucts of datasets from available studies of CP response
%%%All structs have same order of genes for expression (if genes were not
%%%measured assigned value of zero to all samples).
%%%The response (as defined in the relevant study) is used to find positive
%%%(Px) and negative (Nx) samples. 

% [AUCx, scorex, labelx, jj] = classifyImmuneCOMP: is the function used for
% classification based on IMPRES score of checkpoint relations (as done for
% neuroblastoma samples) 

%%%For data from WARGO and TIM CHAN studies, the samples are divided to
%%%datasets based on pre/on treatment, and for WARGO the drug recieved (a-CTLA-4 or a-PD-1) 


