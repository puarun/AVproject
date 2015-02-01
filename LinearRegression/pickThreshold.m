function tArrayAV = pickThreshold()
%This fn. will help you pick the threshold visually using ginput
%   just execute and use cross hair to pick threshold
% Author: Arun.P.U. 8-5-2014
    load m6dB_youden_noEqualDis_O5dB.mat
    for j = 1:24
        plot(Xarray{1,1}{1,j},Yarray{1,1}{1,j})
        hold on
        plot(Xarray{1,2}{1,j},Yarray{1,2}{1,j},'r')
        plot(Xarray{1,2}{1,j},Yarray{1,3}{1,j},'k')
        vline(tArray{1,3}(j));
        title(['Band-',num2str(j)]);
        [x,y]=ginput(1);
        thresh=x+y-1;
        tArrayAV(j)=thresh;
        hold off
    end
end


