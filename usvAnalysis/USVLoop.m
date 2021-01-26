%loop 1 min file 5 times to make 5 mins file
clear all; close all;
files = dir('*.wav');
wavPath =  'C:\Users\Stowers Lab\Desktop\';
totalMice = size(files, 1);
for j = 1:totalMice
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
% [y, Fs]= audioread(files(j).name);
i = audioinfo(waveFile);
totalSec = i.Duration;
startInd = 1;
stopInd = length(y)-1;
y= y(startInd:stopInd);
yLoop = repmat(y,5,1);
filname= [files(j).name];
waveFileLoop = [wavPath, filname, '_5minsloop.wav'];
display(['Saving trimmed WAV file for '])
    audiowrite(waveFileLoop, yLoop, Fs);
end