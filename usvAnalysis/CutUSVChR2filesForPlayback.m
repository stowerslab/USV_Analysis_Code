clear all; close all;
files = dir('*.wav');
wavPath =  'C:\Users\Stowers Lab\Desktop\';
totalMice = size(files, 1);
mouseNums = {'EC8'};
for j = 1:totalMice
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
% [y, Fs]= audioread(files(j).name);
i = audioinfo(waveFile);
totalSec = i.Duration;

startInd = 113*Fs;
stopInd = 173*Fs-1;
yTrim= y(startInd:stopInd);
  display(['Filtering ', waveFile])

    fpass =  40000; % in Hz

    yTrimFilt = highpass(yTrim, fpass, i.SampleRate) ;
yLoop = repmat(yTrimFilt,5,1);
filname= [files(j).name];
waveFileLoop = [wavPath, filname, '_5minsloop.wav'];
display(['Saving trimmed WAV file for '])
    audiowrite(waveFileLoop, yLoop, Fs);
end