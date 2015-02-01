function [ys,CSpline]=reconstructAudio(yn,C1,fs,Tw,Ts)
% Use this program to reconstruct the individual bands and combine them together 
% to create an audio file.
% Author: Arun.P.U
% Date: November 12 2013
   %Define required params to convert vector to frames
    %L = length( yn{1});          % length of the speech signal
    Nw = round( fs*Tw*0.001 );              % frame duration (in samples)
    Ns = round( fs*Ts*0.001 );              % frame shift (in samples)
    M=25;
    %nfft = 2^nextpow2( 2*Nw );              % FFT analysis length
    %K=nfft/2+1;
    %f = linspace( 0, fs*0.5, K ); % frequency range (Hz), size 1xK
    %ys=zeros(1,L);
    ys=0;
    ysub=0;
    % For each subband you can then replicate the mask computed and thereby
    % increase length to FFT resolution. Then multiply with the framed
    % signal.
    ampFac=1;
    [b,c] = designMelfilterbank(M,fs);
    CSpline=zeros(size(C1));
   % C1(C1==0)=0.1;% to vary between 0.5 and 1
    for m=1:M
        x=filtfilt(b(m,:),1,yn);
        %ind=f>=c(m)&f<=c(m+1);
        % divide noisy and clean speech signals into frames
        [ frames.x, indexes ] = vec2frames( x, Nw, Ns, 'rows', @hanning, false );
        %repmat MASKband in order to match length of number of frames
        span = 2;
       % CSpline(:,m)=smooth(C1(:,m),span);% smoothen transitions between 0 and 1
        CSpline(:,m)=C1(:,m);% dont smooth ***
        MASKband=repmat(CSpline(:,m),1,size(frames.x,2));        
        frames.proc = MASKband.*frames.x;
        % perform overlap-and-add synthesis
        y1 = frames2vec( frames.proc, indexes, 'rows', @hanning, 'G&L' ); 
        % truncate extra padding (back to original signal length)
        ys = ys+y1;  
        
%          % reconstruct audio file with subMASK derived from the Oracle
%         subMASK1=repmat(subMASK(:,m),1,size(frames.x,2));
%         frames.sub=ampFac*(frames.x(1:size(subMASK,1),:).*subMASK1);
%         % perform overlap-and-add synthesis
%        % indexes1=indexes(1:size(MASKnew,1),:);
%         y2 = frames2vec( frames.sub, indexes, 'rows', @hanning, 'G&L' ); 
%         % truncate extra padding (back to original signal length)
%         ysub = ysub+y2;        
    end
   %---Visualization--------------------------------------------------%
   supressOutput=0;
   if supressOutput
       % Writing output to files
       %rmsVal=rms(speech.clean);
       %y=(rmsVal/rms(y))*y;
       ys=(rmsVal/rms(ys))*ys;
       %wavwrite(y,fs,'idbmCleanup.wav');
       wavwrite(ys,fs,'subbandCleanup.wav');
       wavwrite(speech.noisy,fs,'noisyAudio')
       figure
       subplot(3,2,1),plot(speech.noisy),title('Noisy speech')
       subplot(3,2,2),myspectrogram(speech.noisy,fs)
       subplot(3,2,3),plot(y),title('IDBM')
       subplot(3,2,4),myspectrogram(y,fs)
       subplot(3,2,5),plot(y),title('Classified')
       subplot(3,2,6),myspectrogram(ys,fs)
   end