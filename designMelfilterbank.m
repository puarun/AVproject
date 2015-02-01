function [b,c]=designMelfilterbank(M,fs)
%This function will desin a mel spaced FIR filter bank. 
%The default sampling frequency is 22050.
%The number of filters, width of transition bands can be set. 
%It uses the fir2 function to do filter design.
%Author: Arun.P.U.
%Date: October 31, 2013

    if ~exist('M','var') M= 25; end% number of filters
    N = 400;%Order of filter, ...
    %limit is number of samples(time_window*nyquist)/3
    % Mel scale warping function
    h2w = @(hz)(1127*log(1+hz/700)); % Hertz to mel warping function
    w2h = @(mel)(700*exp(mel/1127)-700); % mel to Hertz warping function
%%
    %Define the mel scales
    R=[300 3000];
    %fs=44100/2;
    f_min = 0;          % filter coefficients start at this frequency (Hz)
    f_max = 0.5*fs;     % filter coefficients end at this frequency (Hz)
    f_low = 68;       % lower cutoff frequency (Hz) for the filterbank 
    f_high = f_max;      % upper cutoff frequency (Hz) for the filterbank
    f = linspace( f_min, f_max, M ); % frequency range (Hz), size 1xK
%%
    % filter cutoff frequencies (Hz) for all filters, size 1x(M+2)
    %c = w2h( h2w(f_low)+[0:M+1]*((h2w(f_high)-h2w(f_low))/(M+1)) );
    c = w2h( h2w(f_low)+[0:M]*((h2w(f_high)-h2w(f_low))/(M)) );
    c=c/f_max;
%%  
    % Design the filters
    % First filter    
    %w=1/M;
    t=0.5*(c(2)-c(1));
    b=[];
    if M==2
        c=[10 1000 2000]/f_max;
        temp1 = fir1(N,[10/f_max 1000/f_max],'bandpass');
        temp2 = fir1(N,[1000/f_max 2000/f_max],'bandpass');
        b=[temp1; temp2];
        return;
    end  
%     A=[1 1 0 0];
%     F=[c(1) c(2) c(2)+t 1];
%     b1 = fir2(N,F,A);
%     b=[b;b1];
    %fvtool(b1,1)
%%
    % Middle filters
    for i = 1:M-1% use M-1 if f_high=f_max
        if i<4
            t=0.5*(c(i+1)-c(i));
        else
            t=0.1*(c(i+1)-c(i));
        end
        A=[0 0 1 1 0 0];
        F=[0 c(i)-t c(i) c(i+1) c(i+1)+t 1];
        b1 = fir2(N,F,A);
        b=[b;b1];
        %fvtool(b1,1)
    end
    %%
    % Final filter if you are ending at 0.5*fs
    if f_high >= f_max  
        t=0.5*(c(end)-c(end-1));
        A=[0 0 1 1];
        F=[0 c(end-1)-t c(end-1) c(end)];
        b1 = fir2(N,F,A);
        b=[b;b1(1:length(b))];
    end
    % Match magnitude response of all filters
   %standardizeFiltergain()
    %Visualize filter magnitude responses.    
    % for h =1:25
   %      fvtool(b(h,:),1)
   %  end
    % Visualize the filters
    %fvtool(b1,1)
   % fvtool(b1(1:end),1)
    %testMelfilterbank(b,fs,c);
    %testSentence(b,M);
end
function testMelfilterbank(b,fs,c)
    Tw=32e-3;% Analysis window length
    %L=floor(fs*Tw);% Number of samples to be analyzed
    %N=2^nextpow2(L);% Order of filter
    t=0:1/fs:Tw;% Time function
    M=100;
   % y=0;% Initialize
    %f = 2*floor(linspace(20,1e4,10))/fs;
    y=0;
%     for k =2:length(c)
%         y=y+M*sin(pi*(c(k-1)+c(k))*t);
%     end
    k=8+1;% The first digit represents the filter number where 
    %you should see a response
    f=(fs/2)*(c(k-1)+(c(k)-c(k-1))/2)
    y=y+M*sin(2*pi*f*t);
    for i = 1: size(b,1)
        %f=0.5*(c(i)+c(i+1))*fs*0.5;
        figure(i)
        subplot(2,1,1)
        plot(t/fs,y)
        out=filtfilt(b(i,:),1,y);
        title('Sine wave')
        %figure(j)
        subplot(2,1,2)
        plot(t,out)
        w = waitforbuttonpress;
        close all
    end   
end
function testSentence(b,M)
    testfile = '/Volumes/DATAHDD/Past_projects/Portland/AV_project/data/pete rspin wav 2/7-1_84.65.wav';
    [sound,fs]=wavread(testfile);
    sound=decimate(sound,2);
    fs=fs/2;
    rmsVal=rms(sound);
    y=0;
    for i = 1:M
        %filter the audio into subbands
        subbandSpeech=filtfilt(b(i,:),1,sound);
        %combine the audio
        y=y+subbandSpeech;       
    end
    y=(rmsVal/rms(y))*y;
    wavwrite(y,fs,[pwd,'\results\','filterAlone.wav']);
end

    