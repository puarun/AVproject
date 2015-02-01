function visualizeClassifiers(pArray)
%Once the path is entered this routine will visualize all classifiers
%   update path at the appropriate variable.
% Author: Arun.P.U.
    %     path = 'C:\Work\AV_project\code\BinaryMask\idbm\testData\results\sameFPS\';
    %     files=dir([path,'*.mat']);
    %     for i = 1:length(files)
    %         fileName=getfield(files,{i},'name');
    %         snr1=strsplit(fileName,'_');
    %         snr=snr1{1};
    %         filePath=[path,fileName];
    %         load(filePath);
    %         figure(3)
    %         subplot(3,2,i)
    %         plot(AUCArray{1},'*-')
    %         hold on,plot(AUCArray{2},'r*-')
    %         plot(AUCArray{3},'k*-')
    %         title(['AUC-',snr])
    %         ylim([0.5,1]);
    %         xlabel('Frequency band')
    %         ylabel('Percentage')
    %         %plot percorrect
    %         figure(4)
    %         subplot(3,2,i)
    %         plot(pArray{1},'*-')
    %         hold on,plot(pArray{2},'r*-')
    %         plot(pArray{3},'k*-')
    %         title(['Performance-',snr])
    %         %legend('Video','Audio','Audio+Video')
    %         ylim([0.5,1]);
    %         xlabel('Frequency band')
    %         ylabel('Percentage')
    %     end
   figure(1)   
   title('Percent right from cMat')
   %figure(2)
  % plot(AUCArray);
   %title('Area under the curve') 
    for i =1:3
        Noise(:,i)=pArray{i}(:,1);
        Sig(:,i)=pArray{i}(:,2);   
    end    
    subplot(2,1,1)
    plot(Noise,'*-')
    ylim([0,1])
    xlabel('Frequency band')
    ylabel('Noise percent right')
    legend('V','A','A+V')
    subplot(2,1,2)
    plot(Sig,'*-')
    ylim([0,1])
    xlabel('Frequency band')
    ylabel('Signal percent right')
    legend('V','A','A+V')
end
