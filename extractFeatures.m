function [FBE,MASK,ah3] = extractFeatures(speech,Tw,Ts,fs,noise,VisualizeFeatureMASK)
% This function will extract all the 15 features (FBE), MASK 
% and ah3 is just a handle to plot.
% Author: Arun. P.U.

   %% INITIALIZE
   % Explode samples to the range of 16 bit shorts
   % if( max(abs(speech))<=1 ), speech = speech * 2^15; end;
    %speech=speech./max(speech);% Normalize speech
    Nw = round( 1E-3*Tw*fs );    % frame duration (samples)
    Ns = round( 1E-3*Ts*fs );    % frame shift (samples)
    nfft = 2^(nextpow2( Nw )+2);     % length of FFT analysis; we are adding +2 to improve FFT
    LC = 1;% Local SNR criterion use 2 for 0 dB SNR
    % resolution
    % Traingular filter bank params to extract features from the FFT
    % magnitude
    K = nfft/2+1;                % length of the unique part of the FFT
    M =15;                        % number of traingular windows plus one
    R = [0 400];
   % window='hanning';
    %% HANDY INLINE FUNCTION HANDLES
    % Forward and backward mel frequency warping (see Eq. (5.13) on p.76 of [1]) 
    % Note that base 10 is used in [1], while base e is used here and in HTK code
    hz2mel = @( hz )( hz );     % Hertz to mel warping function
    mel2hz = @( mel )( mel); % mel to Hertz warping function
    % Type III DCT matrix routine (see Eq. (5.14) on p.77 of [1])
   % dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
     %                                  .* repmat(pi*([1:M]-0.5)/M,N,1) ) );
    % Cepstral lifter routine (see Eq. (5.12) on p.75 of [1])
    %ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
    %% FEATURE EXTRACTION 
    % Preemphasis filtering (see Eq. (5.1) on p.73 of [1])
    % speech = filter( [1 -alpha], 1, speech ); % fvtool( [1 -alpha], 1 );
    % Framing and windowing (frames as columns)
    %noise=abs(speech-clean);
    framesSignal = vec2frames( speech, Nw, Ns, 'cols', @hamming, false );
    framesNoise = vec2frames( noise, Nw, Ns, 'cols', @hamming, false );
    % Magnitude spectrum computation (as column vectors)
    MagSignal = abs( fft(framesSignal,nfft,1) );
    % Triangular filterbank with uniformly spaced filters 
    H = trifbank( M, K, R, fs, hz2mel, mel2hz ); % size of H is M x K 
    %H(end+1,:)=1;
    % Filterbank application to unique part of the magnitude spectrum
    FBE = H * MagSignal(1:K,:);
  %  FBE=standardizeFeatures(FBE);
    %FBE = H * magnitude;
   % FBE = normalizeFeatures(FBE')';
    %FBE=FBE(1:15,:);
   % FBE=FBE/max(max(FBE));
%________________________________________________________________________________________________________
%% MASK COMPUTATION
%---------------------------------------------------------------------------
%   SNR computation
    %SNR = @(signal,noisy)( 20*log10(norm(signal)/norm(signal-noisy)) );
%   MagNoise = abs( fft(framesNoise,nfft,1) );
    %MagSignal = MagSignal(1:K,:);
    %MagNoise = 
    % magnitude=MAG(1:K,:);
    %magnitude=normalizeFeatures(magnitude);
    % freq = fs/2*linspace(0,1,K);
    % plot(freq(1:100),MAG(1:100,:));
    twoNorm = @(M)(sqrt(sum(abs(M).^2,1))); %# The two-norm of each column
    rmsSignal=twoNorm(framesSignal);
    rmsNoise=twoNorm(framesNoise);
    SNR=20*log10(rmsSignal./rmsNoise); % SNR computation
    MASK=zeros(size(SNR));
    MASK(SNR>LC)=1;
    if VisualizeFeatureMASK
        ah3=subplot(4,1,1);
        plot(SNR)
        ylim([-1 10]);
        hline(LC,'r');
        title('Signal to Noise Ratio')
    else
        ah3=[];
    end
  %  visualizeSignals(speech,noise,SNR);
%     [bl,al]=butter(9,400/fs,'low');
%     envelopeSpeech2=filtfilt(bl,al,envelopeSpeech1);
%     envelopeClean2=filtfilt(bl,al,envelopeClean1);
end
% function FBEnew2=standardizeFeatures(FBE)
%     FBEnew=bsxfun(@rdivide,FBE,FBE(1,:));
%     FBEnew1=FBEnew(2:end,:);
%     mu=mean(FBEnew1,2);
%     sigma=std(FBEnew1')';
%     FBEnew2=bsxfun(@minus,FBEnew1,mu);
%     FBEnew2=bsxfun(@rdivide,FBEnew2,sigma);
%     FBEnew2=FBEnew2';
% end
% function visualizeSignals(speech,noise,SNR)
%    subplot(3,1,1),plot(speech)
%    subplot(3,1,2),plot(noise)
%    subplot(3,1,3),plot(SNR)
%    ylim([0 10]);
%    w=waitforbuttonpress;   
% end
% EOF
