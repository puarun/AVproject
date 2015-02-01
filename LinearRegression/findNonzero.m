function featA = findNonzero(path)
% Scan the beta array and find the non-zero features.
% These non-zero features alone will be used in fitting the non-linear
% model
% Author: Arun.P.U. 8-5-2014
    %load m10dB_pruneAudio.mat
    load(path);
    featA={};
    for j = 1:3
        feat={};
        temp = BeatArray{1,j};
        for band = 1:25            
            [i,j,v]=find(temp(band,:));
            feat{band}=j;
        end
        featA=[featA;feat];
    end
    %save('nonzero_Audio_m10.mat','feat');
end