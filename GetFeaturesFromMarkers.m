function lipFeatures = GetFeaturesFromMarkers( filePath)
% [x, dx, d2x, features] = GetVideoFeaturesFromMarkers( filePath, subtractDrift, relativeTo)
% This function extracts video features that have been recorded by the
% Luxland API.  It includes 66 video features positioned around the mouth,
% chin, eyes, and nose.  This function retrieves the position of these
% markers and arranges them in vectors.  
% x:  This is the position of the marker for different video frames. The
% rows represent the marker and the columns represent the frame.  The first
% 66 rows represent the x-position, and rows 67-132 represent the
% y-position of markers 1-66.
% dx:  The frame-to-frame differences (i.e. derivative of the markers)
% d2x:  The frame-to-frame accelerations 
% features:  These have been extracted based on subjective heuristics -
% primarily that the width and height of the lips, mouth and chin provide
% the most relevant information.  This can later be evaluated by performing
% a max-mutual information analysis on the pixels.  But for now, I am just
% using these heuristics
% tx:  Represents only the x-axis marker values (i.e. 1-66 of x vector)
% ty:  Represents only the y-axis marker values (i.e. 67-132 of x vector)

% We use the foillowing 4 video features: lip height, lip width, standard
% deviation of lip height and standard deviation of lip width.
% Last update:  8-5-2014
% Author: ARUN.P.U.
    %% INITIALIZE
subtractDrift='true';relativeTo=0;fitEllipse=0;
COL_X=1;COL_Y=2;
txtDat = csvread(filePath);
framedata.frames = txtDat(1:length(txtDat(:,1)),1);
num_markers = (length(txtDat(1,:)))/2;
framedata.markers = [txtDat(:,COL_X) txtDat(:,COL_Y)];
for ix=2:num_markers
    framedata.markers = [framedata.markers txtDat(:,COL_X+2*(ix-1)) txtDat(:,COL_Y+2*(ix-1))];
end
%% IF THERE ANY DRIFT SUBRACT HERE
% The eye centers are located at markers 0 and 1.  We will use these to
% remove any offset as the face moves during the recording. 
if subtractDrift==true
    left_eye_start = [mean(framedata.markers(1:3,1)) mean(framedata.markers(1:3,2))];
    right_eye_start = [mean(framedata.markers(1:3,3)) mean(framedata.markers(1:3,4))];
    left_eye_diff = txtDat(:,1:2)-repmat(left_eye_start,length(txtDat(:,1)),1);
    right_eye_diff = txtDat(:,3:4)-repmat(right_eye_start,length(txtDat(:,1)),1);
    mean_diff = [mean([left_eye_diff(:,1) right_eye_diff(:,1)],2) mean([left_eye_diff(:,2) right_eye_diff(:,2)],2)];
    framedata.markers=txtDat-repmat(mean_diff,1,length(txtDat(1,:))/2);
end

x=framedata.markers';
totalMarkers = length(x(:,1))/2;

if relativeTo>0
    relx=x(relativeTo,:);
    rely=x(relativeTo+totalMarkers,:);
    relx=repmat(relx,totalMarkers,1);
    rely=repmat(rely,totalMarkers,1);
    rel=[relx;rely];
    x=x-rel;
end
% Extract the derivative and acceleration of the features
dx(:,1)=zeros(length(x(:,1)),1);
for a=2:length(x(1,:))
    dx(:,a)=x(:,a)-x(:,a-1);
    d2x(:,a)=dx(:,a)-dx(:,a-1);
end

% Extract the heuristic features
tx=x(1:2:132,:);
ty=x(2:2:132,:);
%% CONVERT MARKER DATA TO ARRAYS
% Get feature amplitudes
MOUTH_TOP = 55;
MOUTH_BOTTOM = 56;
MOUTH_LEFT_TOP = 57;
MOUTH_LEFT_BOTTOM = 59;
MOUTH_RIGHT_TOP = 58;
MOUTH_RIGHT_BOTTOM = 60;
MOUTH_TOP_INNER = 62;
MOUTH_BOTTOM_INNER = 65;
MOUTH_LEFT_TOP_INNER = 61;
MOUTH_LEFT_BOTTOM_INNER = 64;
MOUTH_RIGHT_TOP_INNER = 63;
MOUTH_RIGHT_BOTTOM_INNER = 66;
MOUTH_LEFT_CORNER = 5;
MOUTH_RIGHT_CORNER = 4;
CHIN_BOTTOM = 12;
CHIN_LEFT = 10;
CHIN_RIGHT = 11;

features(1,:)=ty(MOUTH_BOTTOM,:)-ty(MOUTH_TOP,:);  
features(2,:)=ty(MOUTH_LEFT_BOTTOM,:)-ty(MOUTH_LEFT_TOP,:);
features(3,:)=ty(MOUTH_RIGHT_BOTTOM,:)-ty(MOUTH_RIGHT_TOP,:);  
%features(4,:)=ty(MOUTH_BOTTOM_INNER,:)-ty(MOUTH_TOP_INNER,:);    
%features(5,:)=ty(MOUTH_LEFT_BOTTOM_INNER,:)-ty(MOUTH_LEFT_TOP_INNER,:);  
%features(6,:)=ty(MOUTH_RIGHT_BOTTOM_INNER,:)-ty(MOUTH_RIGHT_TOP_INNER,:);  
features(4,:)=tx(MOUTH_LEFT_CORNER,:)-tx(MOUTH_RIGHT_CORNER,:);  
features(5,:)=tx(MOUTH_RIGHT_TOP,:)-tx(MOUTH_LEFT_TOP,:);
features(6,:)=tx(MOUTH_RIGHT_BOTTOM,:)-tx(MOUTH_LEFT_BOTTOM,:);

feat(1,:)=cartdist(tx(MOUTH_BOTTOM,:),ty(MOUTH_BOTTOM,:),tx(MOUTH_TOP,:),ty(MOUTH_TOP,:));
feat(2,:)=cartdist(tx(MOUTH_LEFT_BOTTOM,:),ty(MOUTH_LEFT_BOTTOM,:),tx(MOUTH_LEFT_TOP,:),ty(MOUTH_LEFT_TOP,:));
feat(3,:)=cartdist(tx(MOUTH_RIGHT_BOTTOM,:),ty(MOUTH_RIGHT_BOTTOM,:),tx(MOUTH_RIGHT_TOP,:),ty(MOUTH_RIGHT_TOP,:));
feat(4,:)=cartdist(tx(MOUTH_LEFT_CORNER,:),ty(MOUTH_LEFT_CORNER,:),tx(MOUTH_RIGHT_CORNER,:),ty(MOUTH_RIGHT_CORNER,:));
feat(5,:)=cartdist(tx(MOUTH_RIGHT_TOP,:),ty(MOUTH_RIGHT_TOP,:),tx(MOUTH_LEFT_TOP,:),ty(MOUTH_LEFT_TOP,:));
feat(6,:)=cartdist(tx(MOUTH_RIGHT_BOTTOM,:),ty(MOUTH_RIGHT_BOTTOM,:),tx(MOUTH_LEFT_BOTTOM,:),ty(MOUTH_LEFT_BOTTOM,:));

% for i =1:6
%     feat(i,:)=(feat(i,:)-min(feat(i,:)))/(max(feat(i,:))-min(feat(i,:)));
% end
%% DEFINITION OF LIP WIDTH AND HEIGHT BASED ON MEAN VALUES OF LEFT, RIGHT AND CENTER
lipWidth=mean(feat(4:6,:));
lipHeight=mean(feat(1:3,:));
%Area=pi*lipHeight.*lipWidth;
lipFeatures1=[lipHeight;lipWidth];
lipFeatures=addVidfeatures(lipFeatures1); %Add in the standard deviation
%%
%features(8,:)=ty(CHIN_BOTTOM,:)-ty(MOUTH_TOP,:);  
%features(9,:)=ty(CHIN_RIGHT,:)- ty(MOUTH_RIGHT_TOP,:);
%features(10,:)=ty(CHIN_LEFT,:) - ty(MOUTH_LEFT_TOP,:);

% % Get features of an ellipse around the outer mouth
%  xouter(1,:)=tx(MOUTH_TOP,:);youter(1,:)=ty(MOUTH_TOP,:);
%  xouter(2,:)=tx(MOUTH_LEFT_TOP,:);youter(2,:)=ty(MOUTH_LEFT_TOP,:);
%  xouter(3,:)=tx(MOUTH_RIGHT_TOP,:);youter(3,:)=ty(MOUTH_RIGHT_TOP,:);
%  xouter(4,:)=tx(MOUTH_BOTTOM,:);youter(4,:)=ty(MOUTH_BOTTOM,:);
%  xouter(5,:)=tx(MOUTH_RIGHT_BOTTOM,:);youter(5,:)=ty(MOUTH_RIGHT_BOTTOM,:);
%  xouter(6,:)=tx(MOUTH_LEFT_BOTTOM,:);youter(6,:)=ty(MOUTH_LEFT_BOTTOM,:);
%  xouter(7,:)=tx(MOUTH_RIGHT_CORNER,:);youter(7,:)=ty(MOUTH_RIGHT_CORNER,:);
%  xouter(8,:)=tx(MOUTH_LEFT_CORNER,:);youter(8,:)=ty(MOUTH_LEFT_CORNER,:);
% 
% % Get features of an ellipse from the upper outer mouth to the chin
%  xchin(1,:)=tx(MOUTH_TOP,:);ychin(1,:)=ty(MOUTH_TOP,:);
%  xchin(2,:)=tx(MOUTH_RIGHT_CORNER,:);ychin(2,:)=ty(MOUTH_RIGHT_CORNER,:);
%  xchin(3,:)=tx(MOUTH_LEFT_CORNER,:);ychin(3,:)=ty(MOUTH_LEFT_CORNER,:);
%  xchin(4,:)=tx(CHIN_BOTTOM,:);ychin(4,:)=ty(CHIN_BOTTOM,:);
%  xchin(5,:)=tx(CHIN_LEFT,:);ychin(5,:)=ty(CHIN_LEFT,:);
%  xchin(6,:)=tx(CHIN_RIGHT,:);ychin(6,:)=ty(CHIN_RIGHT,:);
%  
% cols=length(features(1,:));
% 
% %h = figure;
% for ix=1:cols
% %      hold all;
% %      scatter(xouter(1,ix),youter(1,ix));hold all;
% %      scatter(xouter(2,ix),youter(2,ix));hold all;
% %      scatter(xouter(3,ix),youter(3,ix));hold all;
% %      scatter(xouter(4,ix),youter(4,ix));hold all;
% %      scatter(xouter(5,ix),youter(5,ix));hold all;
% %      scatter(xouter(6,ix),youter(6,ix));hold all;
% %      scatter(xouter(7,ix),youter(7,ix));hold all;
% %      scatter(xouter(8,ix),youter(8,ix));hold all;     
% %      legend('top','right corner','left corner','chin bottom','chin left','chin right');
%     X=[xouter(:,ix) youter(:,ix)]';
%     minX=min(X(1,:));maxX=max(X(1,:));
%     minY=min(X(2,:));maxY=max(X(2,:));
%     [mouth_ellipse long_axis short_axis] = compute_guaranteedellipse_estimates(X);
%     features(11,ix)=long_axis;
%     features(12,ix)=short_axis;
% %    a = mouth_ellipse(1); b = mouth_ellipse(2); c = mouth_ellipse(3);
% %    d = mouth_ellipse(4); e = mouth_ellipse(5); f = mouth_ellipse(6);
% 
% %    fh = @(x,y) (a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f);
% %    h = ezplot(fh,[minX maxX minY maxY]);
% %    axis([minX maxX minY maxY]);  
% %    pause;
% %    mouth_ellipse=fit_ellipse(xouter(:,ix),-1*youter(:,ix));
% %    chin_ellipse=fit_ellipse(xchin(:,ix),-1*ychin(:,ix));
% %    features(13,ix)=chin_ellipse.short_axis;
% %    features(14,ix)=chin_ellipse.long_axis;
%     
% end

% Get feature derivatives
% cols=length(features(1,:));
% rows=length(features(:,1));
% features((rows+1):rows*2,:)=zeros(rows,cols);
% features((rows+1):rows*2,2:cols)=features(1:rows,2:cols)-features(1:rows,1:(cols-1));
% cols=length(feat(1,:));
% rows=length(feat(:,1));
% feat((rows+1):rows*2,:)=zeros(rows,cols);
% feat((rows+1):rows*2,2:cols)=feat(1:rows,2:cols)-feat(1:rows,1:(cols-1));
%% NORMALIZE THE VIDEO FEATURES
    lipFeatures=lipFeatures';
    temp=bsxfun(@minus,lipFeatures,min(lipFeatures));%minus min
    temp=bsxfun(@rdivide,temp,max(lipFeatures)-min(lipFeatures));% divide by max minus min
    lipFeatures=temp;
    lipFeatures=lipFeatures';
end
function d=cartdist(x1,y1,x2,y2)
% some data
    % x=[0,3,4,5];
    % y=[0,4,5,6];
% the engine
     d=sqrt((x1-x2).^2+(y1-y2).^2);
end
