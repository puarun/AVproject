function C2=testPerformance(featset,BeatArray,tArray,SNR,meaA1,sigA1,filepath,coeff,doPCA,MASKa,MASKv,vName,fps)
   % featset='video';
   % Beta=B45_nostd;   
   %Flags
    if ~exist('writeFile','var');writeFile=1;end
    if ~exist('doPCA','var');doPCA=1;nfu=15;end
    rocCurve=0;
    visualizeSeparation=0;
    %%
    %*********************
    % Get Beta values here
    %**********************
    if length(BeatArray)>1        
        switch featset
           case 'video'
               indx=1;
           case 'audio'
               indx=2;
            case 'audiovideo'
                indx=3;
        end
        Beta=BeatArray{1,indx};
        tA=tArray{1,indx};
        meaA=meaA1{1,indx};
        sigA=sigA1{1,indx};
    else
        Beta=BeatArray{1};
        tA=tArray{1};
        meaA=meaA1{1};
        sigA=sigA1{1};
    end
    M=25;
    %%
    %*****************************
   % All the signal processing gets done here       
    %Get audio and video features as a combined array
    %***************************88
    [vFeatures,audioStart,audioEnd]=videoFeat(vName);
    [FeatureArray1,MASK1,Tw,Ts,fs,c,yn,y,yclean]=extractMASKnfeatures(filepath,...
        [],SNR,vFeatures,M,audioStart,audioEnd,fps);
    switch featset
        case 'audiovideo'
            FeatureArray=FeatureArray1;
        case 'audio'
            FeatureArray=FeatureArray1(:,1:45,:);% for 45 features 45!
        case 'video'
            FeatureArray2=FeatureArray1(:,46:end,:);
            %*************************
            %Additional video features
            addVideo=0;
            if addVideo
                FeatureArray = addVidfeatures(FeatureArray2);
            else
                FeatureArray = FeatureArray2;
            end
            %*************************
    end
    %%
    %*********************
    % Classification gets done here    
    %***********************
    C1=[];
    auc=[];
    auc1=[];
    for j = 1:M-1
        bj=Beta(j,:)';
        FeatTest1=FeatureArray(:,:,j);
        me=meaA(j);
        si=sigA(j);
        FeatTest=conditionFeatures(FeatTest1,me,si);
        if doPCA
            coeff1=coeff{j};
            if j==1
                FeatTest1=coeff1*FeatTest(:,1:30)';
            else
                FeatTest1=coeff1*FeatTest';
            end
            %prune features to use
            p = glmval(bj,FeatTest1(1:nfu,:)','logit');
        else
            p = glmval(bj,FeatTest,'logit','constant','off');
        end            
        C1=[C1,p];
        if visualizeSeparation
            subplot(2,1,1),scatter(p,MASK1(j,:))
            title(['Band - ',num2str(j)])
            subplot(2,1,2),hist(p,100)
            waitforbuttonpress
        end
        if rocCurve                
            yj=MASK1(j,:);
            try 
                [X,Y,T,AUC] = perfcurve(yj,p,1);
%                 plot(X,Y)
%                 title(['Band - ',num2str(j),' AUC - ',num2str(AUC)])
%                 waitforbuttonpress
            catch
                AUC=0;
            end
            auc=[auc,AUC];
        end
    end
    auc1=[auc1;auc];
    %  C1=C1';
    % Make the binayr masks
    C2=makeMASK(C1,tA);
    C2(:,end+1)=zeros;
    %******************************
    % Modifying masks by combining audio and video classifiers
    %****************************
%     if strcmp(featset,'audiovideo')
%         C2=combineMasks(MASKa,MASKv,C2);
%     elseif strcmp(featset,'audio')
%      %   C2=combineMasks(C2,MASKv,[]);
%     end
%%
%****************************************
    % Reconstruct audio signal
%*******************************
    [ys,C21]=reconstructAudio(yn,C2,fs,Tw,Ts);
    % Write output files
    if writeFile            
        %ys=SSBoll79(ys,fs,0.45); % spectral subtraction, removes silent noise from first 0.45 s
        rmsVal=rms(yn);
        ys=(rmsVal/rms(ys))*ys;
        y=(rmsVal/rms(y))*y;
        temp=strsplit(filepath,'\');
        fName=temp{6};        
        switch featset
            case 'video'
                wavwrite(ys,fs,[pwd,'\results\','Video_',num2str(SNR),'_',fName]);                
            case 'audio'
                wavwrite(ys,fs,[pwd,'\results\','Audio_',num2str(SNR),'_',fName]);                
            case 'audiovideo'
                wavwrite(ys,fs,[pwd,'\results\','AudioVideo_',num2str(SNR),'_',fName]);
                wavwrite(yn,fs,[pwd,'\results\','Noisy_',num2str(SNR),'_',fName]);
                wavwrite(y,fs,[pwd,'\results\','IDBM_',num2str(SNR),'_',fName]);
        end           
    end

    %C3=C21';C4=flipud(C3);C5=fliplr(C4);
    %MASK3=flipud(MASK1);%MASK4=fliplr(MASK3);
%     subplot(2,1,2),imshow(~MASK3),title('Ideal Binary MASK')
%     subplot(2,1,1),imshow(~C5),title(['Predicted Binary MASK with ',featset,' features'])
    % Recurse only once    
end
function C2=makeMASK(C1,tA)
    C2=C1;
    for k = 1:size(C1,2)       
        thresh=tA(k);
        C2(C2(:,k)<=thresh,k)=0;
        C2(C2(:,k)>thresh,k)=1;
    end
end
function newFeat=conditionFeatures(TestFeat,mu,sigma)
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
 %addpath('C:\Work\AV_project\code\BinaryMask\idbm');
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