function runGMMBayes(directory,snr)
    % This script is used to train the GMM-Bayes classifier.
    % Before executing this script please make sure all paths are 
    % set correctly. The resulting values are saved into the file
    % at the very end.
    % Author: Arun.P.U.
    %clear all; close all; clc; randn('seed',0); rand('seed',0); fprintf('./n');
    %% load up all the data
    if ~exist('directory','var') || isempty(directory)
        directory='/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression/results/61714/'; % Feature path for training set
    end
    addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression');% Linear model fit path
    %snr=-4;
    load ([directory,'Features-',num2str(snr),'dB-51-feat-nostd']);% load up features
    %load ('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/results/Feature_sets/Test_WhiteNoise_-10dB-51-feat-nostd.mat');
    load ('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/results/Oracles/Oracles_5dB-51-feat-nostd.mat');% load up oracles
    %load nonzero_Audio_m10.mat
    %% Make the test and train features the same  
    TestFeat = TrainFeat;
    totalMASK=[TrainMASK;TestMASK];% using 100% of sentences for training
    TrainMASK = totalMASK;
    TestMASK = totalMASK;
    bayesSA={};
    cmatA={};
    meaAA={};
    sigAA={};
    postprobA={};
    %% Train and test the classifiers
    path = [directory,'Regression-',num2str(snr),'dB_pruneAll.mat']; % Regression model path
    feat = findNonzero(path);% find non-zero features from the linear model
    for i =1:3 % loop for video, audio and audio+video
        feat1=feat(i,:);
        [bayesS,cmat,meaA,sigA,postprob] = featureAnalysisGMMBayes( TrainFeat,...
            TrainMASK,TestFeat,TestMASK,feat1); % call function to fit the model
        bayesSA=[bayesSA;bayesS];% Gaussian Mixture Model - Bayes Classifier
        cmatA{i}=cmat;% confusion matrix
        meaAA{i}=meaA;% mean
        sigAA{i}=sigA;% sigma
        postprobA{i}=postprob;% posterior probability
    end
    %% save the results
    readme=sprintf('bayesSA -  array contaning GMM-Byaes models/ncmatA - confusion matrix array/nmeaAA - mean matrix array/nsigA - standard deviation array/npostprob - posterior probability array/nAll arrays have 3 dimensions, /none each for video, audio and AV in that order.');
    save([directory,'GMMBayes-',num2str(snr),'dB'],'bayesSA','cmatA','meaAA','sigAA','postprobA','readme');

end