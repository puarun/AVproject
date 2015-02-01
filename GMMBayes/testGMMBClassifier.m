function C2=testGMMBClassifier(bayesS,SNR,meaA1,sigA1,filepath,vName,fps,k)
% This function will make use of the generated GMM-Bayes model and produce
% the output files. The input parameters are the following:
% bayesS: GMM-Bayes model
% SNR: signal to noise ratio
% meaA1: mean Array
% sigA1: standard devaitation array
% filpath: the file to which all three model need to be applied
% corresponding video path
% fps : frame rate of video file
% k: video:1; audio:2; audio-video: 3
% C2: output MASK
% Author: Arun.P.U.
   % featset='video';
   % Beta=B45_nostd;   
   %Flags
   addpath('gmmbayestb-vOriginal/');
   addpath('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression');
    if ~exist('writeFile','var');writeFile=1;end
    if ~exist('doPCA','var');doPCA=1;nfu=15;end
    rocCurve=0;
    visualizeSeparation=0;
    M=25;
%% All the signal processing gets done here; get audio and video features as a combined array
    
    vFeatures=GetFeaturesFromMarkers(vName);
    audioStart=[];audioEnd=[];
    [FeatureArray1,MASK1,Tw,Ts,fs,c,yn,y,yclean]=extractMASKnfeatures(filepath,...
        [],SNR,vFeatures,M,audioStart,audioEnd,fps);    
    FeatureArray=FeatureArray1;   
%% Apply the GMMBayes Classifier that was trained   
   % SNR1 = -1*SNR;
  % type='Audio';% audio or AV
    C1=[];
    %load (['nonzero_',type,'_m',num2str(SNR1),'.mat']); % What features to use for which band?
    pathlin = ['/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression/results/White_Noise/m',num2str(SNR),'dB_pruneAll.mat'];%linear classifier path
    feat1 = findNonzero(pathlin);
    feat = feat1(k,:);
    for j = 1:M-1
        feature=feat{j};
        FeatTest1=FeatureArray(:,feature,j);
        me=meaA1{j};
        si=sigA1{j};
        FeatTest=conditionFeatures(FeatTest1,me,si);
        result = testClassifier(FeatTest,bayesS{j});
        C1=[C1,result];
       % disp(['Confusion matrix for band ',num2str(j),' is: ']);
       % disp(num2str(confusionmat(MASK1(j,:)+1,result)));
    end
    %C1=C1';  
    C2=C1;
    C2(:,end+1)=zeros;
    %% Visualize the mask
%     imagesc(flipud(C2'));
%     colorbar();
%     waitforbuttonpress;
%% Reconstruct audio signal
    [ys,C21]=reconstructAudio(yn,C2,fs,Tw,Ts);
    % Write output files
    if writeFile            
        %ys=SSBoll79(ys,fs,0.45); % spectral subtraction, removes silent noise from first 0.45 s
        rmsVal=rms(yn);
        ys=(rmsVal/rms(ys))*ys;
        y=(rmsVal/rms(y))*y;
        temp=strsplit(filepath,'/');
        fName=temp{6};
        switch k
            case 1                
                wavwrite(ys,fs,[pwd,'/results/','Video_',num2str(SNR),'_',fName]);
                wavwrite(yn,fs,[pwd,'/results/','Noisy_',num2str(SNR),'_',fName]);
                wavwrite(y,fs,[pwd,'/results/','IDBM_',num2str(SNR),'_',fName]);
            case 2
                wavwrite(ys,fs,[pwd,'/results/','Audio_',num2str(SNR),'_',fName]);
            case 3
                wavwrite(ys,fs,[pwd,'/results/','AudioVideo_',num2str(SNR),'_',fName]);
        end                
    end

    %C3=C21';C4=flipud(C3);C5=fliplr(C4);
    %MASK3=flipud(MASK1);%MASK4=fliplr(MASK3);
%     subplot(2,1,2),imshow(~MASK3),title('Ideal Binary MASK')
%     subplot(2,1,1),imshow(~C5),title(['Predicted Binary MASK with ',featset,' features'])
    % Recurse only once    
end
%%helper functions
function result = testClassifier(Ptest,bayesS)
    % This is the Bayesian case.
    pdfmat = gmmb_pdf(Ptest, bayesS);% probability density function
    postprob = gmmb_normalize( gmmb_weightprior(pdfmat, bayesS) );% compute posterior probability 
    temp = postprob(:,1);
    adjustprob = 0;% boost signal probability by 0.8
    postprob(:,2) = postprob(:,2)+adjustprob; % Use the probability values for making mask 
    result = gmmb_decide(postprob); % which posterior probability is greater? Signal or noise?
    result(result==0)=1; % making it the regular way *** % this is due to weird 
    %line in the toolbox -Arun P.U.  
  %  disp(['Zeroing ',num2str(total),' elements'])
    % Make the binary masks
    result=result-1; % making mask regular way ***
    probMASK = 0;% apply probability values instead of MASK
    if probMASK        
        zeroindex = result==0;
        result(zeroindex) = temp(zeroindex);
    end
end
function newFeat=conditionFeatures(TestFeat,mu,sigma)
% this function applies the saved mu and sigma to the features of 
% file being processed
    newFeat=[];
    for i = 1:size(TestFeat,3)        
        FBE=TestFeat(:,:,i);
        FBEnew2=bsxfun(@minus,FBE,mu);
        sigma(sigma==0)=1;
        FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
        newFeat=cat(3,newFeat,FBEnew2);
    end
end
function newFeat=standardizeFeatures(TestFeat)
    newFeat=[];
    for i = 1:size(TestFeat,3)        
        FBE=TestFeat(:,:,i);
        mu=mean(FBE);
        sigma=std(FBE);
        FBEnew2=bsxfun(@minus,FBE,mu);
        sigma(sigma==0)=1;
        FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
        newFeat=cat(3,newFeat,FBEnew2);
    end
end
function [vFeatures,audioStart,audioEnd]=videoFeat(vName)
    vFeatures=GetFeaturesFromMarkers(vName);
    %**********************          
    % Trim audio and video
    %***********************
    trimFeatures=0;
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
    %*******************************
end