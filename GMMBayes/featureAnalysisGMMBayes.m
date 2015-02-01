function [bayesS,cmat,meaA,sigA,postprobA] = featureAnalysisGMMBayes( TrainFeat,TrainMASK,TestFeat,TestMASK,feat1 )
% This function will make use of the GMMBayes toolbox to fit a model.
%   We start off by adding a path to the toolbox
    %loadOracleFeatures
    addpath('gmmbayestb-vOriginal\');
meaA={};
sigA={};
cmat={};
P=[];
% switch type
%     case 'video'
%         TrainFeat=TrainFeat1(:,46:end,:);
%         TestFeat=TestFeat1(:,46:end,:);
%     case 'audio'
%         TrainFeat=TrainFeat1(:,1:45,:);
%         TestFeat=TestFeat1(:,1:45,:);        
%     case 'audiovideo'
%         TrainFeat=TrainFeat1;
%         TestFeat=TestFeat1;
% end
%% Fit a Gaussian mixture model for each band with two classes (Noise - 0; Signal - 1)
    for band = 1:24
        if band ==1 && 1
            feat=feat1{1,band};
            %feat=46:51;
        else
            %feat=feat{1,band};
            feat=feat1{1,band};
        end 
        TestFeatBand = TestFeat(:,feat,band);
        TrainFeatBand = TrainFeat(:,feat,band);
        %********************************************
        %Normalize features and save meand and std
        %********************************************
        [mea,sig,TrainFeatBand1]=standardizeFeatures(TrainFeatBand);
        TestFeatBand1=conditionFeatures(TestFeatBand,mea,sig);
        meaA{band}=mea;
        sigA{band}=sig;
        %********************************
        %Define test and train sets
        %***********************************
        Ptrain=TrainFeatBand1;
        Ttrain=TrainMASK(:,band)+1;
        Oraclet=TestMASK(:,band)+1;
        XXt=TestFeatBand1;                   
        %**********************************************
        % Prepare the training data set for each class
        %**********************************************
        disp(['Starting training for band ',num2str(band)])                       
        bayesS{band} = trainClassifier(Ptrain,Ttrain,band,feat);            
        %****************************************
        %Prepare test data set & test classifier           
        %Bayesian classifier; Apply Bayes formula
        %****************************************
        if ~isempty(bayesS{band})            
            disp(['Starting testing for band ',num2str(band)])         
            [OracleP,postprob] = testClassifier(XXt,bayesS{band});
            postprobA{band}=postprob;
            disp(['Confusion matrix for band ',num2str(band),' is: ',])
            cmat1=confusionmat(Oraclet,OracleP)
            cmat{band}=cmat1;
        else
            disp('Empty Bayes sructure ...')
            cmat1=[];
            cmat{band}=cmat1;
        end
    end
%% Plot the confusion matrix and save the GMM just created.
   % plotCmat(cmat);
  %  save('GMMBayes_Audio_-10dB_GEM.mat','meaA','sigA','cmat','bayesS');
end
function bayesS = trainClassifier(data,class,band,feat)
    %*************************
    %Compute kmeans
    %*************************    
    nComp=2*length(feat);
   % nComp=2;
%     N=size(Ptrain,2);
%     if nComp<N
%         nComp=N;
%     end
%% Define all initial parameters
   % FJ_params = { 'Cmin',2,'Cmax', nComp, 'thr', 1e-6, 'animate', 1 ,'verbose', true};
    %EM_params = { 'init', 'fcm1', 'components', nComp, 'thr', 1e-8,'verbose',true};
    GM_params = {'Cmax', nComp, 'animate', true,'verbose',true};    
    try
%% Train the classifier by calling the appropriate function in the toolbox
        disp(['Running GEM with max components set to ',num2str(nComp)]);
        %bayesS = gmmb_create(data, class, 'FJ', FJ_params{:});        
        bayesS = gmmb_create(data, class, 'GEM', GM_params{:});
%         bayesS = gmmb_create(data,class,'EM',EM_params{:} );
%     catch
%         try
%              disp('Second try ..');
%                 % FJ_params{2} = 200;          
%                 % bayesS = gmmb_create(data, class, 'FJ', FJ_params{:});
%                   GM_params{1}= 200;
%                    bayesS = gmmb_create(data, class, 'GEM', GM_params{:});
% %                 EM_params{2}=200;
% %                 bayesS = gmmb_create(data,class,'EM',EM_params{:} );        
%         catch
%             try
%               disp('Third try ..');
%                 % FJ_params{2} = 200;          
%                 % bayesS = gmmb_create(data, class, 'FJ', FJ_params{:});
%                   GM_params{1}= 200;
%                    bayesS = gmmb_create(data, class, 'GEM', GM_params{:});
% %                 EM_params{2}=100;
% %                 bayesS = gmmb_create(data,class,'EM',EM_params{:} );      
%                 
%             catch
%                 try
%                     disp('Fourth try ..');
%                     % FJ_params{2} = 200;          
%                     % bayesS = gmmb_create(data, class, 'FJ', FJ_params{:});
%                      GM_params{1}= 200;
%                       bayesS = gmmb_create(data, class, 'GEM', GM_params{:});
% %                     EM_params{2}=50;
% %                     bayesS = gmmb_create(data,class,'EM',EM_params{:} );      
                catch
                    bayesS=[];
%                 end
%             end
%         end
    end
        
end
function [result,postprob] = testClassifier(Ptest,bayesS)
%% Test the classifier with the test data set
    pdfmat = gmmb_pdf(Ptest, bayesS); % Fit a probability density function.
    postprob = gmmb_normalize( gmmb_weightprior(pdfmat, bayesS) ); % Multiply by the priors
    result = gmmb_decide(postprob);% compare probability of belonging to noise vs signal class
    result(result==0)=1;% Assign class labels
    total=sum(result==0);
    disp(['Zeroing ',num2str(total),' elements'])
end
function covar = computeCovar(cluster1,cluster2,ctrs)
    for i = 1:size(ctrs,2)        
        cluster1(:,i)=cluster1(:,i)-ctrs(1,i);
        cluster2(:,i)=cluster2(:,i)-ctrs(2,i);        
    end
    %covar = (cluster1*cluster2')/(length(cluster1)-1);
    covar1 = cov(cluster1);
    covar2=cov(cluster2);
    covar = cat(3,covar1,covar2);
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
function [mu,sigma,newFeat]=standardizeFeatures(TrainFeat)
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
function plotCmat(cmat)
    bands=size(cmat,2);
    temp1=zeros(1,bands);
    temp2=zeros(1,bands);
    for k = 1:bands
        if ~isempty(cmat{k})
            temp1(k)=cmat{k}(1,1)/(cmat{k}(1,1)+cmat{k}(1,2));
            temp2(k)=cmat{k}(2,2)/(cmat{k}(2,1)+cmat{k}(2,2));
        else
            temp1(k)=0;
            temp2(k)=0;
        end
    end
    plot(temp1,'Marker','o','MarkerFaceColor','b','MarkerSize',10),...
        xlabel('Band'),ylabel('Percent right')
        hold all 
    plot(temp2,'Marker','d','MarkerFaceColor','r','MarkerSize',10),...
        xlabel('Band'),ylabel('Percent right')
    title('Classification Rate of Signal and Noise'),ylim([0,1.1]);
    legend('Noise','Signal');hold off,
end