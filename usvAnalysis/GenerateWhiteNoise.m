clear all; close all;
wavPath =  'C:\Users\Stowers Lab\Desktop\';
noise = wgn(11520000,1,75); % generate 1 mins of white noise @ 75dB
%noise = rand(11520000,1); % can't speicific the power. Need to measure for
%real if want to use rand

Fs=192000;
waveFile = [wavPath, 'whitenoise_60s_75dB.wav'];
    audiowrite(waveFile, noise, Fs);
