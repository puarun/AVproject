function  plotPercorrect(pArray )
% Given the percentage right this function will plot 
% according to each class (noise and signal) the percentage
% of samples identified correctly. This visualization is just a convinient
% way of evaluating the results.
% Author: Arun.P.U. 8-5-2014
temp=[];
    for i = 1:3         
        temp=[temp,pArray{i}(:,1)];        
    end
    subplot(2,1,1)    
    plot(temp,'*-')
    ylim([0,1])
    ylabel('Noise percent right')
    legend('V','A','AV')
    
temp=[];
    for i = 1:3                
        temp=[temp,pArray{i}(:,2)] ;        
    end
    subplot(2,1,2)
    plot(temp,'*-')
    ylim([0,1])
    ylabel('Signal percent right')
     legend('V','A','AV')    
    title('Classification accuracy for White Noise and SNR = -8dB SPL')    
end