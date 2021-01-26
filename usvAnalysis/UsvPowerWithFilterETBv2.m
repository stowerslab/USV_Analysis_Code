% top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior & plot power is USV band
% NOTE: requires Statistics toolbox

clear all; close all;
% wavPath = 'C:\data\Jingyi\Behavior\2018-2-15 ChR2 T2 USV\';
 %wavPath =  'E:\Emily\Behavior data\2018-9 ArCHRT Esr-Cre and Controls\ArchT sorted by light\All light OFF\';
% wavPath = 'E:\Emily\Behavior data\2018-9 ArCHRT Esr-Cre and Controls\ArchT sorted by light\All light ON\';
%wavPath =   'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Female T1 T2 n=8\50Hz\';
%wavPath = '/run/media/emily/BALTZ/'; 
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Terminal stimulation PAG n=5\50hz\';
wavPath = 'C:\Users\Stowers Lab\Desktop\';

files = dir('*.wav');

totalMice = size(files, 1);
   
%mouseNums = {'FEC3_T1_column3', 'FEC3_T1_column10'};
% mouseNums = {'EC10_T1'};
%totalMice = size(mouseNums, 2); 

usvPwrPerSampleAll = [];
timepointUSVs = [];


for j = 1:totalMice
   % mouseNum = mouseNums{1,j};

% mouseNum = 'litter1_uc2'; 
% mouseNum = 'litter1_uc3'; 
% mouseNum = 'litter1_uc4'; 
% mouseNum = 'litter3_uc3'; 
% mouseNum = 'SeaWave _20180926 140133'; 
waveFile = [wavPath, files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
% [y, Fs]= audioread(files(j).name);
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
% note that P is the power spectral density in W/Hz
Tperiod = i.Duration/length(T);

refPower = 10^-12; %reference power is 10?12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 

disp(['Filtering noise...'])
% idea here take z-score & compare to other freqs to remove broadband noise: 
zsignal = zscore(signal);
lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
highFreq = find(F>80000,1,'first'); % index for high cutoff used below
zsignal(1:lowFreq,:) = 0;           % lowpass - set everything below cutoff to zero
zsignal(highFreq:end,:) = 0;        % highpass - add by Jingyi on 11/29/2018
zthresh = 1.5;                      % 1.3 for pup USV
zsignal(zsignal<zthresh) = 0;       % threshold zscore
signalCleaned = signal;             % create a copy to clean below
signalCleaned(zsignal==0) = 0;      % JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

% calculate power in the whistle snippets over time
disp(['Calculating acoustic power...'])
usvPowerPerSample = mean(signalCleaned); % average power across frequencies (so mean dB from ~40-80kHz, over temporal smoothing filter); do NOT normalize for now

% for filtering/smooothing:
wndsz = round(0.05/Tperiod); % convert seconds to samples
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 
AUCpower(j) = trapz(usvPowerPerSample)/totalSec; 
% AUCpowerARRAY(i) = AUCpower; 
% usvPwrPerSampleAll(end+1,1:22400) = usvPowerPerSample(1:22400); %ignore:
% this was written in case sampling is off between subjects
    
% change minpeakdistance depending on your expected duration of USVs
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvs(j) = length(locUsvs); %number of USVs

avgPeaks(j)= mean(pks); %average power of USVs per trial
stdPeaks(j) =std(pks); %standard deviation of USV power per trial

timepointUSVs= T(locUsvs); %timepoints that the USVs occurred at
if locUsvs >1 
    firstcall= timepointUSVs(1); 
else
    firstcall= NaN ; 
end
latencytocall(j) = firstcall-5; %(in seconds) 
interUSVinterval= diff(timepointUSVs); %find the difference between timepoints 
avgIUI(j) = mean(interUSVinterval);  %average the inter USV interval per trial
stdIUI(j) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial
    
     
% disp(['Plotting...'])
% % lowFreq = 1; 
% % highFreq = length(F);
% fig=figure;
% 
% subplot(3,1,1)
% title('original');
% imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
% set(gca,'YDir','normal') %flip to make small freq on bottom
% colormap(1-bone);
% % caxis([thresh thresh+30]);
% % colorbar
% ylabel('frequency (kHz)');
% xlabel('time (seconds)')
% 
% subplot(3,1,2)
% title('filtered, z-scored, and thresholded');
% imagesc(T, F(lowFreq:highFreq)./1000, signalCleaned(lowFreq:highFreq,:)); 
% set(gca,'YDir','normal') %flip to make small freq on bottom
% colormap(1-bone);
% % caxis([zthresh zthresh+5]);
% % colorbar
% ylabel('frequency (kHz)');
% xlabel('time (seconds)')
% 
% subplot(3,1,3)
% title('power per sample');
% % plot(T(1:end-1), dUsvDt)
% plot(T(1:validSamples), usvPowerPerSampleSmooth)
% ylabel('USV power (total dB in 40-80kHz band)');
% xlabel('time (seconds)')
% 
% %add by Jingyi, save fig
% 
% filname= [files(j).name];
% saveas(fig, [wavPath, filname, 'USVPower.tif'], 'tif');

end

AUCAvg= mean(AUCpower); 
% 
% hHeatmap = figure;
% heatmapOld(usvPwrPerSampleAll);
% set(gca,'XTick',[find(T>5,1,'first') find(T>10,1,'first')])
% set(gca,'XTickLabel',{'0' num2str(stimLength)})
% x1 = [find(xDown2>0,1,'first') find(xDown2>0,1,'first')];
% x2 = [find(xDown2>stimLength,1,'first') find(xDown2>stimLength,1,'first')];
% yLims = get(gca, 'YLim');
% y1 = [yLims(1) yLims(2)];
% hold on;
% plot(x1, y1, 'b--', 'LineWidth', 2); %line for stimulus start
% plot(x2, y1, 'b--', 'LineWidth', 2); %line for stimulus endcolormap(1-gray);
% caxis([0 50]);
% xlabel('time (sec)', 'FontSize', fontSz);
% title('USV power', 'FontSize', fontSz);
% colorbar
% hold off;
%     
%     
