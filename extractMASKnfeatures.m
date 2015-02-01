function [FeatureArray1,MASK,Tw,Ts,fs,c,yn,y,yclean1]=extractMASKnfeatures(Noise,SpeechNoisy,SNR,vFeatures,M,aSt,aE,fps)
%This script will generate features from the sound file. 
% In addition it also calls some functions written by Kamil Wojcicki, 2011.
% Author: Arun.P.U. % November 7 2013
    
    %clear all; close all; clc; randn('seed',0); rand('seed',0); fprintf('.\n');
    if strcmp(class(Noise),'char')
        [speech,fs,sentenceName,yclean1]=prepareSpeech(Noise,SNR,aSt,aE);
    else
        speech.noisy=SpeechNoisy;% noisy speech
        speech.noise=Noise;% noise
        fs=11025;
        sentenceName='Test';
    end
    %% INITIALIZE
    % Define window parameters
    %---------------------------------------------------------------------------------%
    fpsm=fps;%82.9;
    Tw = (2/fpsm)*1e3;    % analysis frame duration (ms) 
    Ts = Tw/2;  % analysis frame shift (ms)
    %LC = SNR;    % local SNR criterion (LC)
    if ~exist('M','var') M= 25; end% number of filters % Number of FIR filters with mel spacing
    %For test data use 3 and for regular data use 25
    Nw = round( 1E-3*Tw*fs );    % frame duration (samples)
    Ns = round( 1E-3*Ts*fs );    % frame shift (samples)
    [b,c] = designMelfilterbank(M,fs); % call filterbank to get parameters for filtering
    c=c*fs*0.5; % convert filter bank coeeficients from gradients to Hz using sampling frequency
    VisualizeFeatureMASK=0; % this will visualize the MASK as it is being populated
    includeVideo=1;% this will include the video features while returning the feature array
   %% SUBBAND FILTERING AND FEATURE EXTRACTION
   %-------------------------------------------------------------------------------------% 
   % Filter the audio with the filterbank and extract features
   %Normalize speech signal
    %speech.noisy=speech.noisy/max(speech.noisy);
    %speech.clean=speech.clean/max(speech.clean); 
    features=[];
    MASK=[];
    %featuresdT=[];
    %featuresdK=[];
    %MASK=[];
    y=0;
    yclean=0;
    %aSt1=floor(aSt/4);
    %aE1=floor(aE/4);
    for i = 1:M
        %filter the audio into subbands
        subbandSpeech=filtfilt(b(i,:),1,speech.noisy);
        subbandClean=filtfilt(b(i,:),1,yclean1);
        subbandNoise=filtfilt(b(i,:),1,speech.noise);
        
        %Eliminate edge effects by triiming
        subbandSpeech(1:800)=0;%eliminate edge effect
        subbandNoise(1:800)=0;%eliminate edge effect
        
       % subbandSpeech(end-1500:end)=[];%eliminate edge effect
       % subbandNoise(end-1500:end)=[];%eliminate edge effect
       
     %   subbandSpeech = subbandSpeech1(aSt1:aE1);
      %  subbandNoise = subbandNoise1(aSt1:aE1); 
        
        %subplot(2,1,1),plot(subbandSpeech),subplot(2,1,2),plot(subbandClean)
        %subbandClean{i}=filter(b(i,:),1,speech.clean);
        %extract the envelope, decimate by a factor of 3 
        envelopeSpeech1=envelope(subbandSpeech,fs);
        %envelopeClean1=abs(hilbert(subbandClean));
        envelopeNoise1=envelope(subbandNoise,fs);
        %envelope=decimate(envelope1,3);
        %extract features
        [FBEr,subMASK,ah3] = extractFeatures(envelopeSpeech1,Tw,Ts,fs,envelopeNoise1,VisualizeFeatureMASK);% need to change fs if you decimate
        if rms(subbandSpeech)<rms(speech.noisy)*0.01% if the signal becomes too small 
           subMASK=zeros(size(subMASK));% it ends up becoming comparable to noise thus boosting up the SNR
        end
        FBE=FBEr; % these are features 1-15
        %compute delta time
        featuresdT=diff(FBE,1,2);
        %set first del tau to zero
        append=zeros(size(featuresdT(:,1))); % these are the temporal difference features 16-30        
        featuresdT=[append,featuresdT];
        %append to above feature set
        FBE=[FBE;featuresdT];
        %compute del frequency
        if i==1
            append1=zeros(size(FBEr));
            FBE=[FBE;append1];
        else
            L=size(FBEr,1);
            previousBand=features(:,1:L,i-1)';
            featuresdK =FBEr-previousBand;% these are the frequency difference features 31-45
            %append to above feature set
            FBE=[FBE;featuresdK]; % these are all the features with a dimension of 45 (15+15+15)
            % and the loop will create 45 features for each band
        end
        % append new page for the current subband in the master feature
        % array
        %FBE=standardizeFeatures(FBE);
        FBE=FBE';
        %FBE=FBE(:,1:15);%use only the first 15 features
        features=cat(3,features,FBE); % 3D concatenation of the features
        % total no. of frames (all sentences ~ 300*400) * no. of features
        % (45 audio + 4 video) * number of bands (25)
        % Classify the features to generate the subMASK
        %subM=classifyfeatures(FBE);
        MASK=[MASK;subMASK];% The MASK has a dimension of 
        % total no. of frames (all sentences ~ 300*400) * number of bands (25)
        
        %Reconstruct speech with idbm
        % convert vector to frames
        [y1,indexes1]=vec2frames( subbandSpeech, Nw,Ns, 'cols', @hamming, 'false' );
        subMASK1=smooth(subMASK);
        subMASK1=subMASK1'; % we do a smoothening so that the idbm doesnt have musical noise
        % this step is optional
        multMASK=repmat(subMASK1,size(y1,1),1);
        ysubFrames=multMASK.*y1;
        % convert frames to vector
        ysub=frames2vec(ysubFrames, indexes1, 'rows', @hanning, 'G&L');        
        y=y+ysub;
%% VISUALIZE THE MASK AND FEATURES AS THE ARRAY IS BEING POPULATED        
        if VisualizeFeatureMASK
            ah1=subplot(4,1,2:3);plot(FBE);xlabel('Frame number'),ylabel('Feature Value')
            ylim([0 3]);
            title(['Band ',num2str(i),' Sentence ',sentenceName{6}])
            %use sentenceName{6} in actual situation
            ah2=subplot(4,1,4);spy(subMASK,25)
            linkaxes([ah1 ah2 ah3],'x');
            w=waitforbuttonpress;
        end
    end
%    y=y*(rms(speech.clean)/rms(y));
    %sound(y,fs);
    yn=speech.noisy;
    yclean=yclean+subbandClean;
    %%
    %----------------------------------------------------------------------------------%
    %This routine will convert the video features to the same length as the
    %audio then convert them to frames so that it has the same length as
    %the audio frames.
   %Add in the video features
   %Add the videofeatures to the feature array
    if includeVideo
        L1=size(vFeatures,2);
        L2=size(features,1);
        vFeat=[];
        nvf=size(vFeatures,1);
        for h =1:nvf
            xA=1:L1;
            yA=vFeatures(h,:);
            xB=1:L1/L2:L1;
            if length(xB)<L2,
                diffr=L2-length(xB);
                xB=[xB,linspace(xB(end),xB(end)+diffr,diffr)];
            end%if length lesser than intended
            if length(xB)>L2,
                diffr=length(xB)-L2;
                xB1=xB(1:end-diffr);
                xB=xB1;
            end%if length greater than intended
            if xB(end)>xA(end),xB(end)=xA(end);end
            vFeat=[vFeat;interp1(xA,yA,xB)];
        end
        %% include diff
     %   vFeat=vFeat';
      %  vFeat1=diff(vFeat);
      %  vFeat2=[zeros(1,nvf);vFeat1];
%        vFeat2=standardizeFeatures(vFeat2')';
 %       vFeat=standardizeFeatures(vFeat')';
        %vFeat3=[vFeat,vFeat2];
        %scale to audio features
        %maxAudio=max(max(max(features)));
        %     maxVfeat2=(maxAudio./max(vFeat2));
        %     C2 = bsxfun(@times,maxVfeat2,vFeat2);
        %     maxVfeat=(maxAudio./max(vFeat));
        %     C1 = bsxfun(@times,maxVfeat,vFeat);
       %  vFeat3=[vFeat,vFeat2];
        %     vFeat3=normalizeFeatures(vFeat3');
        %includeVideo=1;%Include video features? 
        %% for all bands
        vFeat3=repmat(vFeat',[1,1,25]);
        FeatureArray1=cat(2,features,vFeat3);
    else
            FeatureArray1=features;
    end    
end
function [speech,fs,sentenceName,clean]=prepareSpeech(file1,SNR,aSt,aE)
% This function will prepare speech with the noise type specified. By
% default it is set to broadband white noise.
    sentenceName=strsplit(file1,'\');
    %file.clean='clean.wav';
    % read clean speech samples from wav file 
    [ clean1, fs1, nbits ] = wavread(file1);
    clean=SSBoll79(clean1,fs1); 
    if ~exist('aSt','var') || isempty(aSt) aSt=1; end;
    if ~exist('aE','var') || isempty(aE) aE=length(clean); end;    
    clean=clean(aSt:aE);
    clean = decimate(clean,4);    
    % read noise samples from wav file 
    %[ noise, fs, nbits ] = wavread( 'ssn.wav' );
    
    %decimate audio
    % Specify the noise type here:
    noiseType=2;
    switch noiseType
        case 1
            noiseFile='/Volumes/DATAHDD/Past_projects/Portland/CDA/data/NST Noise.wav';
            [ noise, fs, nbits ] = wavread( noiseFile );
        case 2
            noiseFile='/Volumes/DATAHDD/Past_projects/Portland/CDA/data/White Noise.wav';
            [ noise, fs, nbits ] = wavread( noiseFile );
        case 3
            %*********DONT USE******************
            noiseFile='/Volumes/DATAHDD/Past_projects/Portland/CDA/data/White Noise.wav';
            [ noise, fs, nbits ] = wavread( noiseFile );
            noise=filterWNoise(noise);            
        case 4
            disp('Not adding any noise');
            fs=11025;            
        otherwise
            disp('Error');
    end
    if noiseType ~= 4        
        noise = decimate(noise,2);
        fs=fs/2;
        %clean=SSBoll79(clean,fs); 
        %clean1=clipSpeech(clean);
        clean=clean/norm(clean);
        [speech.noisy,speech.noise] = addnoise( clean, noise, SNR );
    else
        speech.noisy=clean;
        speech.noise=clean*0.001;
    end    
    %noise=filterWNoise(noise1);
    %fs=fs/2;    
    %noise=decimate(noise,2);
    %SNR = 10; % noisy speech signal-to-noise ratio (dB)   
    % mix clean speech and noise at a prescribed SNR    
    %speech.noisy = filter([1 -0.95],1,speech.noisy);
    %speech.clean = filter([1 -0.95],1,speech.clean);
   % rmsVal=rms(speech.noisy);
    % create time vector
  %  time = [ 0:length(clean)-1 ]/fs;
end

function noise2=filterWNoise(noise)
% if you wish to use narrow band noise then use this function
    n = 100; Wn = 0.1;
    ftype = 'low';
    % Transfer Function design
    [b,a] = fir1(n,Wn,ftype);
    noise1=filtfilt(b,a,noise);
    noise2=noise1(300:end);
    
end
function env=envelope(signal,fs)
% this will compute the envelope of the speech signal
    type='hilbert';
    switch type
        case 'hilbert'
           env=abs(hilbert(signal));
        case 'lowpass'
           [bl,al]=butter(5,400/fs);
           env=filtfilt(bl,al,signal);
        case 'env_secant'
            env=env_secant();
    end           
end
%EOF----------------------------------------------------------------------------------------%
