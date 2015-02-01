function makeFeatures(Snratio)
% This function will make features at the specified SNR (Snratio) and 
% noise type (white noise by default). The path to save the files is specified at the end of the function.
% The .mat file saved will have two variables TrainFeat and TestFeat. The
% cutting point is specified by percent (P) of data for training.
% Implemented after Kim et al. JASA 126(3): pp 1486.
% Author: Arun.P.U.

%clear all; close all; clc; randn('seed',0); rand('seed',0); fprintf('./n');
%% Note on how to generate the Oracle
% You can generate the oracle by splitting TrainMASK into TrainMASK and TestMASK
% similar to how we do it for the features useing the crunchFeatset function.
% Date: November 11, 2013
    %% Define paths and what sentences to train on
    %percent of files on which the training is to be done
    P=90;
    %SNR to be trained/tested on
    SNR=Snratio;
    % Specify folder path where you have the sound files
    path='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin wav 2';
    % Read all wave files
    files=dir([path,'/*.wav']);
    s=size(files);
    %N=2;
    %From the percentage input compute the number of files to be read
    %N=0.01*P*s(1);
    % The random file indices to be read is stored in r
    %r=floor(s(2)+(s(1)-s(2)).*rand(N,1));
    %r=1:400;% This specifies which files will train the classifier. 
    %To do all files, set r=1:400 - Arun P.U.
    r=1:10;
    % Initialize MASK and feature array
    TrainFeat=[];
    TrainMASK=[];
    %Video files extraction
    %Include video features
    videopath='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin markers 2';
   % vidFiles=dir([videopath,'*.csv']);
    %% Make feature array and MASK (Oracle)
    % Loop through all files to be read and extract the features and
    % corresponding ORACLE which is a MASKnew (25*No. of frames) derived from
    % the IDBM(nFFT*No. of frames).
    trimFeatures=0;
    for i = 1:length(r)
        %Get audio files
        fName=getfield( files,{r(i)},'name');
        filepath=[path,'/',fName];
        %Get video features
        VF=strsplit(fName,'_');
        vName=[videopath,'/',VF{1},'.csv'];
        t=strsplit(VF{2},'.');
        fps=str2num([t{1},'.',t{2}]);
        vFeatures= GetFeaturesFromMarkers(vName);        
        % Trim audio and video
         if trimFeatures            
            [audio,fs]=wavread(filepath);
            [audioStart,audioEnd,videoStart,videoEnd]=vadSpeechfile(audio,fps,fs);
            if videoEnd>length(vFeatures)
                 videoEnd=length(vFeatures);
            end
            vFeatures=vFeatures(:,videoStart:videoEnd);
         else
             audioStart='';audioEnd='';videoStart=[];videoEnd=[];             
         end        
        %Extract MASK and audio features
        [features,MASKnew,dump1,dump2,fs,dump3,dump4,y,yclean]=extractMASKnfeatures(filepath,[],SNR,vFeatures,25,audioStart,audioEnd,fps);
        y=y*10;
%***************************************************        
%         sound(y,fs)%sound check
%         %visual check
%         figure(1),subplot(2,2,1),plot(y),subplot(2,2,3),plot(yclean)
%         subplot(2,2,2),myspectrogram(y,fs)%,[18,1],@hamming, 1024, [-45 -2], false, 'default', false, 'lp');
%************Visualize mask and spectrogram
%          subplot(2,1,1),myspectrogram(yclean,fs)
%          ylim([0,5500]);
%          M1=flipud(MASKnew);
%          subplot(2,1,2),imshow(M1);
%          waitforbuttonpress
%****************************************************        
        TrainFeat=[TrainFeat;features];% put all features together
        % Make sure MASK and features agree on number of frames
        NframesMASK=length(MASKnew);
        NframesFeatures=size(features,1);
%         if NframesFeatures<NframesMASK
%             Masknew1=MASKnew(1:NframesFeatures,:);
%         else
         Masknew1=MASKnew';
%         end
        % For all the sentences we keep appending MASKnew that contains 25 frequency bands
        % to this array. 
        TrainMASK=[TrainMASK;Masknew1];
        %Frame video signals and make sure the frame length of video and
        %audio agree
        if r(i)==P*0.01*length(r)
            len=length(TrainMASK);
        end
        
    end
    % Start timing the classification
  % tic;
    %%
%     if sum(MASK(:,end))==0
%         count=24;
%     else
%    count=25;
%     end
   % addVideo=1;
  %  TrainFeat=standardizeFeatures(TrainFeat);
    %[O1,dump]=TestTrain('train',TrainFeat,TrainMASK,[]);
    %% divide the whole set into train and test feat
    [TrainFeat,TestFeat]=crunchFeatset(TrainFeat,len); % If you want to traing and test with all 400 sentences dont call
    [TrainMASK,TestMASK]=crunchMASK(TrainMASK,len);
    % this function.
    readme = sprintf('Dimension of features: no. of frames * number of features * number of bands/n Dimension of features: no. of frames * no. of bands/nTestFeat: Test feature set/n TrainFeat: Training feature set/nTrainMASK: Traning MASK');
    save(['Features-',num2str(SNR),'dB','-51-feat-nostd','.mat'],'TrainFeat','TestFeat','TrainMASK','TestMASK','readme');
   %% just checking
  % save(['Features-',num2str(SNR),'dB','-51-feat-nostd','.mat'],'TrainFeat','TrainMASK');
    %cMat1=testScript(O1,SNR);
    %TimeS=toc;
    %load NSTNoise_0dB_AV.mat;
    %save NSTNoise_0dB_AV.mat
    % end timing
end
% function cMat1=testScript(O1,SNR)
%   % Call the test script
%     %t=1:400;
%     %p=setdiff(t,r);
%     p=281:400;
%     %p=282:283;
%     cMat1=TestClassifer(O1,p,SNR);
% end
% function newFeat=standardizeFeatures(TestFeat)
%     newFeat=[];
%     for i = 1:size(TestFeat,3)        
%         FBE=TestFeat(:,:,i);
%         mu=mean(FBE);
%         sigma=std(FBE);
%         FBEnew2=bsxfun(@minus,FBE,mu);
%         sigma(sigma==0)=1;
%         FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
%         newFeat=cat(3,newFeat,FBEnew2);
%     end
% end
function [TrainFeat,TestFeat]=crunchFeatset(TrainFeat,len)
% this function will cut the features into two arrays (TrainFeat and TestFeat)
% at the specified point which is len.
    %load Features-0dB-51-feat-nostd.mat
    TestFeat=TrainFeat(len+1:end,:,:);
    TrainFeat1=TrainFeat(1:len,:,:);
    clear TrainFeat
    TrainFeat=TrainFeat1;
    %save('Features-0dB-51-feat-nostd.mat','TestFeat','TrainFeat');
end
function [TrainMASK,TestMASK]=crunchMASK(TrainMASK,len)
% this function will cut the features into two arrays (TrainFeat and TestFeat)
% at the specified point which is len.
    %load Features-0dB-51-feat-nostd.mat
    TestMASK=TrainMASK(len+1:end,:);
    TrainMASK1=TrainMASK(1:len,:);
    clear TrainFeat
    TrainMASK=TrainMASK1;
    %save('Features-0dB-51-feat-nostd.mat','TestFeat','TrainFeat');
end
%EOF--------------