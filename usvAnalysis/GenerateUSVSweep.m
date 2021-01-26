clear all; close all;
wavPath =  'C:\Users\Stowers Lab\Desktop\';
Fs=192000; % match the ChR2 and WT playback
t = 0:1/Fs:3; % 3s of USV Sweep
Sweep = chirp(t,17000,3,20000); %From 17kHz to 20kHz
Sweep = Sweep';
% pspectrum(Sweep,Fs,'spectrogram','TimeResolution',1, ...
%     'OverlapPercent',99,'Leakage',0.85)
yLoop = repmat(Sweep,20,1);

 waveFile = [wavPath, 'USVsweep.wav'];
audiowrite(waveFile, yLoop, Fs);