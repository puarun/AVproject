function [X1,Y1,AUC1,Beta,tArray,Parray,meaA,sigA,T1,coeffA,yHat]=checkPerformance(featset,TrainFeat,TestFeat,TrainMASK,TestMASK)
% This function will fit a linear model and prune the features 
% if there is a bad fit according to the criteria p<0.05 one at a time.
% Pruning means zeroing out that feature by assigning the corresponding
% beta = 0.
% X1, Y1 - arrays contaning x,y coordinates of the AUc curves
% AUC1 - Area under the cuve for all the bands
% tArray - Threshold picked by the Youden's index
% Parray - percent right correctly identified of the 2 classes (noise (0) and signal (1))
% meaA - mean array
% sigA - standard deviation
% T1 -  threshold array
% coeffA - not relevant if not applying PCA 
% yHat - predicted yHat values for the test feature set
% Author: Arun.P.U. ; 8-5-2014

%% Intialize everything
    %clear all;warning off;    
    %featset='audio';
    doPCA=0;
    addVideo=0;
    nof=size(TrainFeat,2);
    pruneFeat=1;
    switch featset
        case 'video'
            indices=[46:nof];
            if addVideo
                nf=length(indices)+2;% use 51 for all features; 
             %   6 video features+ std height+
            %std width+ entropy height+ entropy width
            % As of August 5 2014 we are using just the first 4 video features
            else
                nf = length(indices);
            end
        case 'audiovideo'
            indices=1:nof;
            nf=length(indices);
        case 'audio'
            indices=1:45;
            % pca
            if doPCA
                nf=15;
            else
                nf=length(indices);
            end
    end
    % initialize beta with dimentions number of bands(25) * number of
    % features (nof)
    Beta=zeros(size(TrainFeat,3),nf);% use nf+1 for a constant
    coeffA=[];
    setEqualDist=1;
    if size(TrainMASK,1)==size(TestMASK,1)
        yHat=zeros(size(TrainMASK));% yhat return values
    else
        yHat=zeros(size(TestMASK));% yhat return values
    end
    for n = 1:size(TrainFeat,3)-1 % we dont train on the last band as it has too 
        % signal and will lead to inconsistent estimate of beta
        yTrain=TrainMASK(:,n);
        yTest=TestMASK(:,n);
        xTrain=TrainFeat(1:length(yTrain),:,n);
        xTest=TestFeat(:,:,n);
        %**********************************
        % In order to test if we trained the signal and noise class with
        % equal number of samples we use the helper function equaldist which will 
        %truncate the additional noise class feature vectors from the training data set 
        if setEqualDist
            [x,y]=equalDist(xTrain,yTrain,n);
            %[xt,yt]=equalDist(xTest,yTest,n);
            xt=xTest;
            yt=yTest;
        else
            x=xTrain;
            y=yTrain;
            xt=xTest;
            yt=yTest;
        end        
        % define train and test set        
        %xtest=x(end-testin:end);
        %ytest=y(end-testin:end);
        %xtrain=x(testin+1:end-testin-1,:);
        %ytrain=y(testin+1:end-testin-1,:);
        %**********************************
        switch featset
            case 'video'
                x1u1=x(:,indices);% the second index is the feature set you want to use
                x2u1=xt(:,indices);
              %  pruneFeat=1;                
                %% Additional video features
                if addVideo               
                   x1u = addVidfeatures(x1u1);                   
                   x2u = addVidfeatures(x2u1);                  
                else
                    x1u=x1u1;
                    x2u=x2u1;
                end
                    
                %% Normalize features
                [mea,sig,xtrain]=standardizeFeatures(x1u);
                xtest=conditionFeatures(x2u,mea,sig);
                ytrain=y;
                ytest=yt;
                %prune features to use
                [X,Y,AUC,beta,thresh,pC,T,yHat1]=pruneFeatures(xtrain,xtest,ytrain,ytest,n,nf,featset,pruneFeat);
            case 'audio'
                x1u=x(:,indices);% the second index is the feature set you want to use
                x2u=xt(:,indices);
                %pruneFeat=1;
                % Normalize features
                [mea,sig,xtrain]=standardizeFeatures(x1u);
                xtest=conditionFeatures(x2u,mea,sig);
                ytrain=y;
                ytest=yt;
                %% Principal Component Analysis
                if doPCA
                    [coeff,xtrain1] = pcaFeatures(xtrain,n);
                    coeffA{n}=coeff;
                    if n==1
                        xtest1=coeff*xtest(:,1:30)';
                    else
                        xtest1=coeff*xtest';
                    end
                    %prune features to use
                    [X,Y,AUC,beta,thresh,pC,T,yHat1]=pruneFeatures(xtrain1(:,1:nf),xtest1(1:nf,:)',...
                        ytrain,ytest,n,nf,featset,pruneFeat);
                else
                    [X,Y,AUC,beta,thresh,pC,T,yHat1]=pruneFeatures(xtrain,xtest,...
                        ytrain,ytest,n,nf,featset,pruneFeat);
                end
           case 'audiovideo'
               pruneSep=0;
                if pruneSep==1
                   %***********************
                    % prune video first
                    %*****************
                    nf=4;
                    featset='video';
                    x1u=x(:,46:45+nf);% the second index is the feature set you want to use
                    x2u=xt(:,46:45+nf);
                   % pruneFeat=1;                
                    % Normalize features
                    [mea,sig,xtrain]=standardizeFeatures(x1u);
                    xtest=conditionFeatures(x2u,mea,sig);
                    ytrain=y;
                    ytest=yt;
                    %prune features to use
                    [X,Y,AUC,betaV,thresh,pC,T,yHat1]=pruneFeatures(xtrain,xtest,ytrain,ytest,n,nf,featset,pruneFeat);
                    %*************************************
                    % prune audio next
                    %***************************
                    nf=45;
                    featset='audio';
                    x1u=x(:,1:45);% the second index is the feature set you want to use
                    x2u=xt(:,1:45);
                   % pruneFeat=1;
                    % Normalize features
                    [mea,sig,xtrain]=standardizeFeatures(x1u);
                    xtest=conditionFeatures(x2u,mea,sig);
                    ytrain=y;
                    ytest=yt;
                    %prune features to use
                    [X,Y,AUC,betaA,thresh,pC,T,yHat1]=pruneFeatures(xtrain,xtest,ytrain,ytest,n,nf,featset,pruneFeat);
                    %********************
                    % put them together
                    %*******************
                    nf = 45+4;
                    featset='audiovideo';
                    indA=find(betaA(1:end-1)~=0);
                    indB=find(betaV(1:end-1)~=0);
                    indB=indB+45;
                    index=1:nf;
                    mea=zeros(1,nf);
                    sig=zeros(1,nf);
                    indices=index([indA,indB]);%use apostrophe if no pruning
                   % nf=length(indices);
                    x1u=x(:,indices);% the second index is the feature set you want to use
                    x2u=xt(:,indices);
                  %  pruneFeat=1;
                else
                    %nf=51;
                    x1u=x(:,1:nf);
                    x2u=xt(:,1:nf);
                   % pruneFeat=1;
                end
                % Normalize features
                [mea1,sig1,xtrain]=standardizeFeatures(x1u);
                xtest=conditionFeatures(x2u,mea1,sig1);
                ytrain=y;
                ytest=yt;
                for g=1:length(indices)
                    h=indices(g);
                    mea(h)=mea1(g);
                    sig(h)=sig1(g);
                end
                %prune features to use
                [X,Y,AUC,beta1,thresh,pC,T,yHat1]=pruneFeatures(xtrain,xtest,ytrain,ytest,n,length(indices),featset,pruneFeat);
                beta=zeros(1,nf);% use 52 for constant
                for j=1:length(indices)
                    val=indices(j);
                    beta(val)=beta1(j);
                end
                %beta(end)=beta1(end);
        end
       
        %############################################################
        %[y,index]=sort(y1);
        %x=x1(index,:);
        %testin=floor(0.05*length(x));
        % define train and test set        
        %xtest=[x(1:testin,:); x(end-testin:end,:)];
        %ytest=[y(1:testin,:); y(end-testin:end,:)];
        %xtrain=x(testin+1:end-testin-1,:);
        %ytrain=y(testin+1:end-testin-1,:);
%         if  sum(ytrain)<sum(ytest)
%             testin_new=floor(sum(ytest)/2);
%             disp(['Not enough ones in this band so using only 2* ',num2str(testin_new),' samples!']);
%             xtest=[x(1:testin_new,:); x(end-testin_new:end,:)];
%             ytest=[y(1:testin_new,:); y(end-testin_new:end,:)];
%             xtrain=x(testin_new+1:end-testin_new-1,:);
%             ytrain=y(testin_new+1:end-testin_new-1,:);            
%         end
        %####################################################################3
        %if AUC>0.5
        X1{n}=X;
        Y1{n}=Y;
        T1{n}=T;
        meaA(n,:)=mea;
        sigA(n,:)=sig;        
        %end
        Parray(n,:)=pC;
        tArray(n)=thresh;
        AUC1(n)=AUC;
        yHat(:,n)=yHat1;
        if n==1
            switch featset
                case 'video'
                     beta1=[beta];% Video only
                case 'audiovideo'
                     beta1=beta;% All 51 features
                case 'audio'
                    if length(beta) ~= nf
                        % remove apostrophe if pruneFeat=1
                        beta1=[beta(1:30),zeros(1,15)];% Audio only us ,beta(31) for const
                    else
                        beta1=beta;
                    end
                    %put ' on beta if pruning 
            end           
            Beta(1,:)=beta1;
        else
            Beta(n,:)=beta;
        end
%         try
%             X1(n,:)=X;
%             Y1(n,:)=Y;
%             AUC(n)=AUC1;
%         catch err
%             X1(n,:)=zeros;
%             Y1(n,:)=zeros;
%             AUC(n)=AUC1;
%         end
%        w=waitforbuttonpress;
    end
%    for i = 1:24,plot(X1(i,:),Y1(i,:)),title([num2str(i),'-',num2str(AUC(i))]),waitforbuttonpress,end
end
function [X,Y,AUC,beta,thresh,pCorrect,T,p]=pruneFeatures(xtrain,xtest,ytrain,ytest,n,nf,featset,pruneFeat)
%% Prune the features based on p < 0.05
%     k = isnan(stats.p);
%     if all(k)
%         return
%     end
    %nf=51;
   % pruneFeat=0;
    xtest1=xtest;
    xtrain1=xtrain;
    ind=[];
    if sum(any(xtest))<nf
        switch featset
            case 'audio'
                ind=1:sum(any(xtest)); % Audio - 45, !30 if using all 45!
            case 'video'
                ind=1:nf; % Video only - 46:51
            case 'audiovideo'
                ind=1:nf;
        end        
    else
        ind=1:nf;
    end
    disp([num2str(size(ind)),' non zero features in band ', num2str(n)])
    xtest1=xtest(:,ind);
    xtrain1=xtrain(:,ind);
    %% intial fit
    [b,dump,stats] = glmfit(xtrain1,ytrain,'binomial','link','logit','const','off');
    % logistic regression
    p = glmval(b,xtest1,'logit','const','off');
    %figure,scatter(y,p);
    %beta=b;
    %beta(k)=0;
    %eliminate the unwanted feature vector
    nF=size(xtest,2);
    %vector=ones(size(xtest,2),1);
%     array=stats.p(1:end-1);
%     nanC=isnan(array);
%     t=find(nanC==1);
    turnoff=[];   
    features=zeros(size(xtest,2),1);  
    features(ind)=1;  
    count=0;
%% pruning stage    
    if pruneFeat
        beta=zeros(1,length(b));
        while max(stats.p)>0.05       
            %trim feature
            array=stats.p;
            [m,ind1]=max(array);            
            xtest1(:,ind1)=[];
            xtrain1(:,ind1)=[];
            ind(ind1)=[];
            % fit the model and calculate p values for fit
            [b,dump,stats] = glmfit(xtrain1,ytrain,'binomial','link','logit','const','off');
            % logistic regression
            p = glmval(b,xtest1,'logit','const','off');        
            %figure,scatter(y,p);        
        end
        disp(['using ',num2str(length(ind)),' features']);
        for j=1:length(ind)
            val=ind(j);
            beta(val)=b(j);
        end
        beta(end)=b(end);
    else
      %  p = glmval(b,xtest1,'logit','const','off'); 
    end
 %% evaulate model performance   
   % beta=b;
   % define class labels for perfcurve
    classLabels=cell(length(ytest),1);
    classLabels(ytest==0)={'Noise'};
    classLabels(ytest==1)={'Signal'};
    % fit probabilities for scores
    [X,Y,T,AUC,OPTR] = perfcurve(classLabels,p,'Signal');
    disp(['band ',num2str(n),' AUC ',num2str(AUC),'-',featset]);
    %indext=find(Y==OPTR(2),1);
    %thresh=T(indext);
    [indext,pCorrect]=pickThreshold(X,Y,T,p,ytest,n);
    thresh=T(indext);    
    if pruneFeat==0
        beta=b;
    end
end
function [mu,sigma,newFeat]=standardizeFeatures(TrainFeat)
% this function computes the mu and sigma from the traning data set
    newFeat=[];
    for i = 1:size(TrainFeat,3)        
        FBE=TrainFeat(:,:,i);
        mu=mean(FBE);
        sigma=std(FBE);
        FBEnew2=bsxfun(@minus,FBE,mu);
        sigma(sigma==0)=1;
        FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
        newFeat=cat(3,newFeat,FBEnew2);
    end
end
function newFeat=conditionFeatures(TestFeat,mu,sigma)
% this function applies the mu and sigma computed from the trainning data
% set to the test data set
    newFeat=[];
    for i = 1:size(TestFeat,3)        
        FBE=TestFeat(:,:,i);
        FBEnew2=bsxfun(@minus,FBE,mu);
        sigma(sigma==0)=1;
        FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
        newFeat=cat(3,newFeat,FBEnew2);
    end
end
function testAllcombinations(xtrain,xtest,ytrain,ytest)
    %initialize optimization parameter
    count=1;
    AUC1=[];
    len=size(x,2);
    for i = 1:3
        temp=nchoosek(1:len,i);
        for j = 1:length(temp) 
            index=temp(j,:);            
            [b,dump,stats] = glmfit(xtrain(:,index),ytrain,'normal');  
            % logistic regression
            p = glmval(b,xtest(:,index),'logit');
            %figure,scatter(y,p);
            % fit probabilities for scores
            [X,Y,T,AUC] = perfcurve(classLabels,p,'Signal');
            %performance(p,y);
           % plot(X,Y)
           % xlabel('False positive rate'); ylabel('True positive rate')
           % title(['ROC for classification by logistic regression'])
      %      w=waitforbuttonpress;
            AUC1=[AUC1,AUC];
        end
        temp1=length(AUC1);
        disp(['Size of array is ',num2str(temp1)]);
        AUCArray{count}=AUC1;
        count=count+1;
        temp2=length(AUCArray);
        disp(['Number of cells is ',num2str(temp2)]);
        AUC1=[];
    end
    %         AUCArray=AUCArray';
%         imagesc(AUCArray)
%         set(gca,'YTick',1:14)
%         set(gca,'XTick',1:15)
%         colorbar
%         grid on
end
function [x,y]=equalDist(xTrain,yTrain,n)
        %equal distribution of ones and zeros in the training set
        [y1,index]=sort(yTrain);
        testin=2*sum(y1);
        y=y1(end-testin:end);
        x1=xTrain(index,:);
        x=x1(end-testin:end,:);
        disp(['Using ',num2str(testin),' train samples for band ',num2str(n)]);        
        %equal distribution of ones and zeros in the test set
%         [y1t,indext]=sort(yTest);
%         testint=2*sum(y1t);
%         yt=y1t(end-testint:end);
%         x1t=xTest(indext,:);
%         xt=x1t(end-testint:end,:);
%         disp(['Using ',num2str(testint),' test samples for band ',num2str(n)]);
end
function [indext,pCorrect]=pickThreshold(X,Y,T,p,ytest,n)
% the default method of picking the threhsold is by computing the Youden's Index
% however you could also do that by picking it by hand from each AUC curve. This function will allow
% you to do that.
    count=0;% Youden's index
    %thType=2;% K-Means
    notSatisfied=1;
    while notSatisfied
        if count ==0
            [maximum,indext]=max(Y-X);
            pCorrect=visualThreshold(X,Y,T,p,ytest,indext,n);
            notSatisfied=0;% uncomment if no visual threshold
        else
            result=input('Are you satisfied (y/n): ','s');
            if strcmp(result,'y')
                notSatisfied=0;
            else
                disp('Pick a threshold on the ROC curve: ')
                [x,y]=ginput(1);
                [temp,indext]=min(abs(Y-y));
                pCorrect=visualThreshold(X,Y,T,p,ytest,indext,n);                
            end            
        end                
        count = count+1;              
    end    
end
function [wcoeff,score,latent,tsquared,explained]=pcaFeatures(xA,i)
%PCA analysis of the features
%   load up the feature und los gehts!
%     directory='C:\Work\AV_project\code\BinaryMask\idbm\results\DiffFPS\';
%     load ([directory,'Features--6dB-51-feat-nostd']);
%     x=TestFeat;
%     xt=TrainFeat;
%     xA=x(:,1:45,:);    
%     [mea,sig,xtrain]=standardizeFeatures(xA);
%     w=[];s=[];l=[];t=[];e=[];    
    if i==1
        [wcoeff,score,latent,tsquared,explained] = pca(xA(:,1:30));

    else
        [wcoeff,score,latent,tsquared,explained] = pca(xA); 
    end    
end
