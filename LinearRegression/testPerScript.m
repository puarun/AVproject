% Use this script to make audio files with the fitted linear model. The
% indices of files to be used is given by r. Configure path before
% execution.
% Author: Arun.P.U.
clear all; close all; clc; randn('seed',0); rand('seed',0); fprintf('.\n');
%load 'C:\Work\AV_project\code\BinaryMask\idbm\testData\results\justAV\m10dB_youden_EqualDis_O5dB_MM.mat';
SNR=0;
%load (['m',num2str(SNR),'dB_pruneAll.mat']);
load (['/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression/results/White_Noise/','m',num2str(SNR),'dB_pruneAll.mat']);
% MM means modified binary mask so audiovideo is not saved
    % Prepare for writing
    if ~exist('results','dir')
        mkdir('results');
    end
    delete('./results/*.wav')    
    %%
    %**********************************
    % Sentence text extraction
    %***************************************
    sentencePath='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/sentences';
    sentenceFiles=dir([sentencePath,'\*.txt']);
    for j = 1:length(sentenceFiles)
        sName=getfield( sentenceFiles,{j},'name');
        stringNum=strsplit(sName,'.txt');
        stringNum1=strsplit(stringNum{1},'-');
        sNo(j,1) = str2num(stringNum1{2});
        sNo(j,2) = j;    
    end
    sNo1=sortrows(sNo,1);
    r=sNo1(1:(end/2)+1,2);
    r=141;
    %%
    %*************************************
    % Loop through all the files
    %******************************************
    addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm');
    path='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin wav 2';
    videopath='/Volumes/DATAHDD/Past_projects/Portland/pete rspin markers 2';
    % Read all wave files
    %SNR=0;
    files=dir([path,'\*.wav']);
    s=size(files);
   % r=391;% trained on indices 1 to 360 so test 360 to 400!
    % Prepare for writing
%     if writeFile && ~exist('results','dir')
%         mkdir('results');
%     end
%     delete('./results/*.wav')       
    for i = 1:length(r)
        fName=getfield( files,{r(i)},'name');
        filepath=[path,'\',fName];
        %Get video features
        VF=strsplit(fName,'_');
        t=strsplit(VF{2},'.');
        fps=str2num([t{1},'.',t{2}]);          
        vName=[videopath,'\',VF{1},'.csv'];        
        MASKv=testPerformance('video',BeatArray,tArray,SNR,mA,sA,filepath,[],0,[],[],vName,fps);%video
        MASKa=testPerformance('audio',BeatArray,tArray,SNR,mA,sA,filepath,[],0,[],MASKv,vName,fps);% audio
        MASKav=testPerformance('audiovideo',BeatArray,tArray,SNR,mA,sA,filepath,[],0,MASKa,MASKv,vName,fps);% audiovideo        
       % MASKav=testPerformanceAV(SNR,r,1,[],[]);% generate mask based on product features
        %MASKav=MASKa & MASKv;        
    end
    %%
    %*************************************
    % Visualization of ourput files
    %*************************************
    visualizeSound=1;
    if visualizeSound    
        audiopath=[pwd,'\results'];
        files=dir([audiopath,'\*.wav']);
        figure(1)
         for i = 1:length(files)
            %Get audio files
            fName=getfield( files,{i},'name');
            [audio,fs]=wavread([audiopath,'\',fName]);
            ax(i)=subplot(length(files),1,i);myspectrogram(audio,fs);%[12 6],...
               % @hamming, 1024, [-45 -2], false, 'default', false, 'lp');% plot spectrogram
            ylim([0, 5000])
           % xlim([0, 3])
            title(fName);
            linkaxes(ax,'x')
         end
    %     figure(2)
    %     fName=getfield( files,{1},'name');
    %     [audio,fs]=wavread([audiopath,'\',fName]);
    %     subplot(4,1,1),myspectrogram(audio,fs);% plot spectrogram
    %     subplot(4,1,2),imshow(~idbm),title('Ideal Binary MASK');
    %     subplot(4,1,3),imshow(~MASKa),title('Predicted Binary MASK with audio features');
    %     %subplot(5,1,4),imshow(~MASKv),title(['Predicted Binary MASK with video features']);
    %     subplot(4,1,4),imshow(~MASKav),title('Predicted Binary MASK with audio and video features');

    end
%%