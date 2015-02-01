function  applyGMMBayes()
% This function will apply the GMMBayes model generated using the toolbox
% to the sound files in order to produce the processed files.
% Kindly ensure all paths are set correctly before execution.
% This function will call testGMMBayesClassifier()
% Author: Arun.P.U.

    %% Apply generated model to sound files
    SNR=-8;% make sure SNR is consistent with below
    load(['/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression/results/White_Noise/m-8-GMMBayes.mat']);% model path 
    % load the saved means, std and classifier
    % Prepare for writing
    if ~exist('results','dir')
        mkdir('results');
    end
    delete('./results/*.wav')    
    
    %% Get the text for sentences 
    sentencePath='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/sentences';% sentence path
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
    r=16:18; % this is index of the sentence that will be used for testing the model
    %% Loop through all the sentence wav files to be tested    
    addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm');% features path
    path='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin wav 2';% audio path
    videopath='/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin markers 2';% video path
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
for k = 1:3    
    bayesS = bayesSA(k,:);
    meaA = meaAA{1,k};
    sigA = sigAA{1,k};        
    for i = 1:length(r)
        fName=getfield( files,{r(i)},'name');
        filepath=[path,'\',fName];
        %Get video features
        VF=strsplit(fName,'_');
        t=strsplit(VF{2},'.');
        fps=str2num([t{1},'.',t{2}]);          
        vName=[videopath,'\',VF{1},'.csv'];        
        MASKv=testGMMBClassifier(bayesS,SNR,meaA,sigA,filepath,vName,fps,k);%video
        %MASKa=testPerformance('audio',BeatArray,tArray,SNR,mA,sA,filepath,[],0,[],MASKv,vName,fps);% audio
        %MASKav=testPerformance('audiovideo',BeatArray,tArray,SNR,mA,sA,filepath,[],0,MASKa,MASKv,vName,fps);% audiovideo        
       % MASKav=testPerformanceAV(SNR,r,1,[],[]);% generate mask based on product features
        %MASKav=MASKa & MASKv;        
    end    
end
    %% Visualize the cleaned up files  
    visualizeSound=0;
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
end

