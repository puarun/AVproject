% This script is to create a demo of the Binary MASK approach to Audio-Visual 
% multiband voice activity detection should give you a basic idea of how the toolbox works and all associated
% function will be called during the demo at some point.
% This script will do the following:
% Demo – August 8, 2014
% Create a results directory where you want to put all files. 
% Make sure to include the SNR and noise type in the name.
% Step1: Make features
% Step2: Fit linear model
% Step3: Process files with linear model
% Step4: Fit GMM-Bayes
% Step5: Process files with GMMBayes
% Author: Arun.P.U.

%% I) Setup new directory to put all files generated in and make the features
snr=-8;
name = 'Demo_White_Noise_-8dB';
if ~exist(name,'dir'),mkdir(name);end % create a new directory to put all results in
% make sure to include the nosie type and SNR in the name.
direc = [pwd,'/',name,'/']; % here is where we will be storing all the output we produce.
makeFeatures(snr); % create features
temp = dir('Features*.mat');
temp1=getfield(temp,'name');
movefile(temp1,[name,'/',temp1])% move the feature mat file to the directory we just created

%% II) Train the linear model
% add a path to the linear model
addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression');
% fit the linear model
BetaArray =   runRegression(direc,[],snr);
%% III) Process files with the linear model we fit above
% This step is not implementd in the demo but if you would like to do
% please refer to testperScript.m under LinearRegression.
%% IV) Train the GMM-Bayes model
% path to the GMM-Bayes model
%addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/GMMBayes');
%runGMMBayes(direc,snr);% Fit the GMM-Bayes model and save it to direc

%% V) Process files with the GMM-Bayes model we fit above
%applyGMMBayes(direc);% apply the GMMBayes model you just generated and save
% the resulting sound files to direc

%% VI) You can explore the TestGUI directory and check out the GUI in order to see how to
% present the audio files.
% All the best! - Arun. P.U.