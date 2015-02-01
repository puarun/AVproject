function plotCmat(cMatA)
% This function will compute the percentage of noise class identified
% correctly and percentage of thye signal class identified correctly for
% all three classifiers (3 curves) across frequency bands.
% Author: Arun.P.U.
    funnoise = @(cmatband)cmatband(1,1)/sum(cmatband(1,:));% inline function for % noise correct
    funsignal = @(cmatband)cmatband(2,2)/sum(cmatband(2,:));% inline function for % signal correct
    for i = 1:3 % loop through A,V, and AV
        cmatband1=cMatA(:,i);
        for j = 1:length(cmatband1{1})            
            pNoise(i,j) = funnoise(cmatband1{1}{:,j});
            pSignal(i,j) = funsignal(cmatband1{1}{:,j});
        end
    end
    
    % Noise Class plot
    subplot(2,1,1)    
    plot(pNoise','*-')
    ylim([0,1])
    ylabel('Noise percent right')
    legend('V','A','AV')
    
    %Signal class plot    
    subplot(2,1,2)
    plot(pSignal','*-')
    ylim([0,1])
    ylabel('Signal percent right')
    legend('V','A','AV')
end