% top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior & plot power is USV band
% NOTE: requires Statistics toolbox

clear all; close all;
wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-9 GCaMP BNST\2018-9-17 GCAMP BNST_T4_RIGHT\';

mouseNum = 'BG6_T4R_FP'; 

waveFile = [wavPath, mouseNum, '.wav'];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
i = audioinfo(waveFile);
totalSec = i.Duration;

% trim file if necessary by % of total length or seconds
% startSec = 5.94;
% stopSec = 6.02;
% startPercent = startSec/totalSec;
% stopPercent = stopSec/totalSec;
% startInd = round(length(y)*startPercent);
% stopInd = round(length(y)*stopPercent);
startInd = 1;
stopInd = length(y);
yTrim = y(startInd:stopInd);

disp(['Taking FFT...'])
nfft = 512;
window = 512;
noverlap = window*0.5;
% thresh = -90; % threshold in decibels;
[~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
%note that P is the power spectral density in W/Hz
Tperiod =i.Duration/length(T);

signal = 10*log10(abs(P)); %convert to dB for acoustic convention (now signal is dB/Hz)

disp(['Filtering noise...'])
% idea here take z-score & compare to other freqs to remove broadband noise: 
zsignal = zscore(signal);
lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
highFreq = find(F>90000,1,'first'); % index for high cutoff used below
zsignal(1:lowFreq,:) = 0; % lowpass - set everything below cutoff to zero
zthresh = 1.5;
zsignal(zsignal<zthresh) = 0; %threshold zscore
% signal(zsignal==0) = 0; %JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

% calculate power in the whistle snippets over time
disp(['Calculating acoustic power...'])
usvPowerPerSample = sum(zsignal); %take sum scross frequencies; do NOT normalize for now
% for filtering/smooothing:
wndsz = round(0.2/Tperiod); %convert seconds to samples % change from 0.05 to 0.2 Jingyi 10/30 to smooth the results
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize
validSamples = length(usvPowerPerSample)-wndsz+1;
 usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 

disp(['Plotting...'])
% lowFreq = 1; 
% highFreq = length(F);
figure;

subplot(3,1,1)
title('original');
imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
set(gca,'YDir','normal') %flip to make small freq on bottom
colormap(1-bone);
% caxis([thresh thresh+30]);
% colorbar
ylabel('frequency (kHz)');
xlabel('time (seconds)')

subplot(3,1,2)
title('filtered, z-scored, and thresholded');
imagesc(T, F(lowFreq:highFreq)./1000, zsignal(lowFreq:highFreq,:)); 
set(gca,'YDir','normal') %flip to make small freq on bottom
colormap(1-bone);
% caxis([zthresh zthresh+5]);
% colorbar
ylabel('frequency (kHz)');
xlabel('time (seconds)')

subplot(3,1,3)
title('power per sample');
% plot(T(1:end-1), dUsvDt)
plot(T(1:validSamples), usvPowerPerSampleSmooth)
ylabel('huh?');
xlabel('time (seconds)')













