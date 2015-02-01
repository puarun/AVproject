function BeatArray =   runRegression(directory,featsetArray,SNratio)
% This functionm will call checkPerformance and regress on TrainFeat 
% BeatArray - Beta coefficients after fitting
% Trainfeat - training data set
% Testfeat - test data set
% featsetArray - type of model (A,V or AV)% leave empty if you want to fit
% all models
% savefile - flag to save the file or not
% Author: Arun.P.U. 8-5-2014
%% INITIALIZE 
useAll = 0;
load ('/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/results/Oracles/Oracles_5dB-51-feat-nostd.mat');
if ~exist('directory','var') || isempty(directory)
    directory='/Volumes/DATAHDD/Past_projects/Portland/AV_project/code/BinaryMask/idbm/LinearRegression/results/White_Noise/';% use difffps for the original features    
    load ([directory,'Features-',num2str(SNratio),'dB-51-feat-nostd']);
    useAll=1;
else
    %TrainMASK = TrainMASK(1:length(TrainFeat),:);
    temp = dir([directory,'Features*.mat']);
    temp1=getfield(temp,'name');
    load ([directory,temp1]);
end

%% WHAT MODELS TO FIT? 
if ~exist('featsetArray','var') || isempty(featsetArray)
    featsetArray{1}='video';
    featsetArray{2}='audio';
    featsetArray{3}='audiovideo';
end
    %% Make the test and train features the same if using 100% to test and train 
    
if useAll 
    TestFeat = TrainFeat;
    totalMASK=[TrainMASK;TestMASK];% use if you are fitting
    TrainMASK = totalMASK;
    TestMASK = totalMASK;
end
     
   %   TestMASK = TrainMASK;% if you are testing
   
    %% Fit models
    for l = 1:length(featsetArray)
        featset=featsetArray{l};
        [X,Y,AUC,Beta,tA,pC,m,s,T,coeffA,yHat1]=checkPerformance(featset,TrainFeat,TestFeat,TrainMASK,TestMASK);
        %[X,Y,AUC,Beta,tA,pC,m,s,T,coeffA]=checkPerformanceAV();
        Xarray{l}=X;
        Yarray{l}=Y;
        AUCArray{l}=AUC;
        BeatArray{l}=Beta;
        tArray{l}=tA;
        pArray{l}=pC;
        mA{l}=m;
        sA{l}=s;
        yHat{l}=yHat1;
    end
    %% plot results
    plotALL=0;
    plot(pC,'*-'),ylim([0,1])
    if plotALL
        %Plot AUC
        figure(7)
        plot(AUCArray{1},'*-')
        hold on,plot(AUCArray{2},'r*-')
       % plot(AUCArray{3},'k*-')
        title('AUC')
        ylim([0,1]);
        xlabel('Frequency band')
        ylabel('Percentage')
        %plot percorrect 0 class
        figure(1)
        plot(pArray{1}(:,1),'*-')
        hold on,plot(pArray{2}(:,1),'r*-')
      %  plot(pArray{3}(:,1),'k*-')
        title('Performance 0 class')
        legend('Video','Audio')
        ylim([0,1]);
        xlabel('Frequency band')
        ylabel('Percentage')
        %plot percorrect 1 class
        figure(9)
        plot(pArray{1}(:,2),'*-')
        hold on,plot(pArray{2}(:,2),'r*-')
     %   plot(pArray{3}(:,2),'k*-')
        title('Performance 1 class')
        legend('Video','Audio')
        ylim([0,1]);
        xlabel('Frequency band')
        ylabel('Percentage')
         %plot combined percentage
        figure(10)
        plot(sum(pArray{1},2),'*-')
        hold on,plot(sum(pArray{2},2),'r*-')
      %  plot(sum(pArray{3},2),'k*-')
        title('Combined performance')
        legend('Video','Audio')
        ylim([1,1.5]);
        xlabel('Frequency band')
        ylabel('Combined Percentage')
    end
%% save results
    savefile = 1;
    if savefile
        readme = sprintf('Dimension of BeatArray: number of features * number of bands\n BeatArray: Regression coefficients after model fitting\ntArray: Threhsold Array \nmA: mean Array \nsA: standard deviation array\ncoeffA: PCA coefficients\n yHat: yhat values\nAUCArray: Area under the curve array');
        save([directory,'Regression-',num2str(SNratio),'dB_pruneAll.mat'], 'BeatArray','tArray',...
            'pArray','mA','sA','coeffA','yHat','AUCArray','readme');    
    end
    visualizeClassifiers(pArray);
end


    