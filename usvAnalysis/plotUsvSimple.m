% top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior
% NOTE: requires Statistics toolbox

clear all; close all;
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_07_JingyiUsv\chr2\';
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_07_JingyiUsv\WT_femaleUrine\';
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2016_07_BnstEsrChR2_USVs\';
wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\CHR2 experiment\2019-6To7 GAQ Vgat Chr2 Data\calling animals_GAQ2 GAQ4 GAQ6\CNO day3\50HZ20S\';
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_07_JingyiUsv\chr2\cut\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Terminal stimulation Amb n=5\T7 dose curve\25hz 10s\';
% wavPath =  'C:\Users\Stowers Lab\Desktop\';
% mouseNum = 'ChR2snippetsCombined60sec'; 
mouseNum = 'GAQ4_C3_50HZ20S_column2'; 
waveFile = [wavPath, mouseNum, '.wav'];
display(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
i = audioinfo(waveFile);
totalSec = i.Duration;

% trim file if necessary by % of total length or seconds
% startSec = 5.94;
% stopSec = 6.02;
% startSec = 328;
% stopSec = 338;
% startPercent = startSec/totalSec;
% stopPercent = stopSec/totalSec;
% startInd = round(length(y)*startPercent);
% stopInd = round(length(y)*stopPercent);
startInd = 1;
stopInd = length(y);
yTrim = y(startInd:stopInd);

display(['Taking FFT...'])
nfft = 512;
window = 512;
noverlap = window*0.5;
% thresh = -90; % threshold in decibels;
[~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh); 

signal = 10*log10(abs(P));

display(['Plotting...'])
% lowFreq = find(F>200,1,'first'); %indeces correspond to ~35-85kHz
% highFreq = find(F>3000,1,'first');
lowFreq = 1; 
highFreq = length(F);
figure;
imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
set(gca,'YDir','normal') %flip to make small freq on bottom
colormap(1-bone);
% caxis([thresh thresh+30]);
ylabel('frequency (kHz)');
xlabel('time (seconds)')

