function x1 = addVidfeatures( x )
%   This function generates additional video features
%   standard deviation of i_2:i+2
%   Since we are using only 4 video features
%   we will use std of height and width
% Author: Arun P.U. ; August 5 2014
    height = x(1,:);
    width = x(2,:);
    height1=zeros(size(height));
    width1=zeros(size(width));
    wide=2;
% standard deviation of height and width    
    for i = wide+1:length(height)-wide
        height1(i) = var(height(i-wide:i+wide));
        width1(i) = var(width(i-wide:i+wide));
    end
% smooth and compute entropy of width
 %   nF=50;% takes plus minus 2 i.e. 5 frames
%    hE=computeEntropy(nF,height);    
% height
  %  wE=computeEntropy(nF,width);    
% prepare array for output by adding new features    
   x1 = [x;height1;width1];
end
% function heightEntropy = computeEntropy(nF,height)
% % This will compute the entropy across specified number of frames
%     vidFeatw2=standardizeArray(height);
%     %vFw3 = smooth(vidFeatw2);
%     heightEntropy=zeros(size(height));
%     for i =  nF+1:length(vidFeatw2)-nF
%         heightEntropy(i) = wentropy(vidFeatw2(i-nF:i+nF),'logenergy');
%     end
% end
% function x1 = makeOutarray(height1,x)
% % This will put the features together with the original array
%     dim=length(size(x));
%     switch dim        
%         case 2
%             x1=cat(2,x,height1);        
%         case 3            
%             temp11=repmat(height1,1,25);  
%             temp21=reshape(temp11,size(height1,1),1,size(temp11,2));
%             x1=cat(2,x,temp21);
%     end    
% end

