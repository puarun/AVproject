function calibration()
% Use this function to generate a sine wave and connect the output to a
% sound meter. That way you can calibrate the master volume to see what
% percent of the master volume corresponds to thye desired SNR.
% Author: Arun.P.U.

    Fs = 44100;        % Samples per second. 48000 is also a good choice
    toneFreq = 1000;  % Tone frequency, in Hertz. must be less than .5 * Fs.
    nSeconds = 2;      % Duration of the sound
    y = sin(linspace(0, nSeconds*toneFreq*2*pi, round(nSeconds*Fs)));
    y=0.5*(y/rms(y));
    sound(y,Fs); % Play sound at sampling rate Fs