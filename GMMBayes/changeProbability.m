% This script was made to check if boosting the posterior probability of
% the signal class led to an improvement in intelligibility.
% Author: Arun.P.U.

clear all;close all;
load(['C:\Work\AV_project\code\BinaryMask\idbm\LinearRegression\results\White_Noise\m-WN--8-GMMBayes-Cmax.mat']);
% load Oracles
load ('C:\Work\AV_project\code\BinaryMask\idbm\results\Oracles\Oracles_5dB-51-feat-nostd.mat');
totalMASK = [TrainMASK;TestMASK];
% add path to the GMM-Bayes toolbox
addpath('gmmbayestb-vOriginal\');
postprob1=postprobA;
% plot before figure
figure(1)
plotCmat(cmatA)
% loop through the different bands and then change the postprob
adjust = 0.6;
for i = 1:3
    for j = 1:24
        Oraclet = totalMASK(:,j)+1;
        %postprob1{1,i}{1,j}(:,1) = postprob1{1,i}{1,j}(:,1);
        postprob1{1,i}{1,j}(:,2) = postprob1{1,i}{1,j}(:,2)+adjust;
        Ptest = postprob1{1,i}{1,j};            
        result = gmmb_decide(Ptest);% compare probability of belonging to noise vs signal class
        result(result==0)=1;% Assign class labels
        %belonging to noise vs signal class
        cmat1=confusionmat(Oraclet,result);
        cmat{1,j}=cmat1;               
    end
    cmatA1{1,i} = cmat;
end
% plot after figure
figure(2)
plotCmat(cmatA1)

