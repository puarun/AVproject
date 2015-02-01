function nocontIndex = separateFiles( sentDir )
% This function will return indices to files that are fully comprehensible 
% and have no context. These indices to files that are fully comprehensible 
% were determined from a listening test taken by John Condon during
% March 2014. These were sentences not degraded by SSBOLL79 function under
% idbm directory (total: len(restdata) = 140 of a possible 200). We had to run our 
% files through the SSBOLL79 function as the initial recordings had a 
% lot of background noise. After this step this function will also separate
% the files into context and no context sentences.
% Author: Arun.P.U.
    nocontIndex=[];
    load restdata.mat % use only those sentences that were not degraded by the processing.
    for i = 1:length(sentDir)
        sentName=getfield( sentDir,{i},'name');
        type1 = strsplit(sentName,'-');
        type2 = strsplit(type1{2},'.');
        type = str2num(type2{1});
        if ~type && sum(find(i==restdata))
            nocontIndex = [nocontIndex;i];
        end
    end
end

