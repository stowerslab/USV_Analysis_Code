clear all; close all;
files = dir('*.wav');
wavPath =  '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior data\USV playback files\For Shawn MoSeq playback\';
totalMice = size(files, 1);
mouseNums = {'EC20_F'};
Fswrite = 192000;
for j = 1:totalMice
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
% [y, Fs]= audioread(files(j).name);
i = audioinfo(waveFile);
totalSec = i.Duration;

startInd = 1*Fs;
stopInd = 600*Fs-1;% cut for 10 mins
yTrim= y(startInd:stopInd);
  display(['Filtering ', waveFile])

    fpass =  40000; % in Hz

    yTrimFilt = highpass(yTrim, fpass, i.SampleRate) ;
% yLoop = repmat(yTrimFilt,5,1); % if wants to loop for the files
filname= [files(j).name];
waveFile= [wavPath, filname, '_cut10mins.wav'];
display(['Saving trimmed WAV file for '])
    audiowrite(waveFile, yTrimFilt, Fswrite);
end