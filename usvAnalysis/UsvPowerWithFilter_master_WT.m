% top-level script to process USV audio data 
% NOTE: requires Statistics toolbox

% change matlab directory to wav file dir, make sure to change
% usvPowerPerSampleSmoothMat to length of wav files as well
clear all; close all;
%wavPath =   'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Female T1 T2 n=8\50Hz\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\BNST ChR2 USV Master dataset\Terminal stimulation PAG n=5\50hz\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2019-2 ChR2 for slice recording animals\2019-2-8 slice recordeing ChR2_T4\25hz 1s\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2019-4 GAD1-3 WAV file only\T1T2\UNILATERAL\50hz Uni\'; 
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2019-4 GAD1-3 WAV file only\SALINE DAY 4\25hz\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2019-4 GAD1-3 WAV file only\CNO DAY 5\25hz\';
% wavPath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\BNST ChR2 USV Master dataset\10Hz_merged\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2019-4 GAD1-3 WAV file only\All Saline 25HZ\';
wavPath = '\\172.29.164.29\jingyi chen\Dropbox\Jingyi Data backup\Voseq data Master folder\WT mice data\WT USV MALE URINE\';
% wavPath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\BNST ChR2 USV Master dataset\10Hz_merged\'


files = dir('*.wav');
plotandsaveindividual = 0; % only if want to plot USVs and save individual Mat files
totalMice = size(files, 1);
%save all meta data for USV analysis  
stat1 = 'AUCpower';  %Not normalized
stat2 = 'NumUSV';
stat3 = 'AvgPeaks';
stat4 = 'stdPeaks';
stat5 = 'latencytoUSV'; %
stat6 = 'AvgIUI'; %
stat7 = 'stdIUI'; % 
stat8 = 'numUsvs0to5'; % 
stat9 = 'numUsvs5to10'; %   
stat10 = 'numUsvs10to15'; % 
stat11 = 'numUsvs15to20'; % 


USVStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[],stat8,[],stat9,[],stat10,[],stat11,[]);

usvPwrPerSampleAll = [];
timepointUSVs = [];
nameList = {};
% usvPowerPerSampleSmoothMat = zeros(1,179963); % for 4 mins WT SMUF analysis only
usvPowerPerSampleSmoothMat = zeros(1,89963); % for 2 mins WT SMUF analysis only
% usvPowerPerSampleSmoothMat = zeros(1,22463);  % hardcode with amount of validSamples %% PETE - 22463 for 30s video
% usvPowerPerSampleSmoothMat = zeros(1,37463);  % hardcode with amount of validSamples % for 50s long video, added by Jingyi on 2019/8/7


for j = 1:totalMice
   % mouseNum = mouseNums{1,j};
waveFile = [wavPath, files(j).name];
nameList = [nameList files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);
% [y, Fs]= audioread(files(j).name);

 % error checking for different sample rates -ZG
    disp(['Checking the sampling rate: ', num2str(Fs)])
    if Fs == 250000
        y = resample(y, 192000, 250000);
    end
    Fs = 192000;
    
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
yTrim = y(startInd:stopInd);%This is for pre-cut file only

disp(['Taking FFT...'])
nfft = 512;
window = 512;
noverlap = window*0.5;
% thresh = -90; % threshold in decibels;
[~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
%note that P is the power spectral density in W/Hz
Tperiod =i.Duration/length(T);

refPower = 10^-12; %reference power is 10?12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 

disp(['Filtering noise...'])
% idea here take z-score & compare to other freqs to remove broadband noise: 
zsignal = zscore(signal);
lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
highFreq = find(F>90000,1,'first'); % index for high cutoff used below
zsignal(1:lowFreq,:) = 0; % lowpass - set everything below cutoff to zero
zsignal(highFreq:end,:) = 0; % highpass - add by Jingyi on 11/29/2018
zthresh = 1.5; %1.3 for pup USV
zsignal(zsignal<zthresh) = 0; %threshold zscore
signalCleaned = signal; %create a copy to clean below
signalCleaned(zsignal==0) = 0; %JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

% calculate power in the whistle snippets over time
disp(['Calculating acoustic power...'])
% usvPowerPerSample = sum(signalCleaned,1); % Use the sum of 40-90kHz energy as USV power per sample
usvPowerPerSample = sum(signalCleaned,1) ./ (highFreq-lowFreq);%average power across frequencies (so mean dB/Hz from ~40-90kHz, over temporal smoothing filter); do NOT normalize for now
% usvPowerPerSample = mean(signalCleaned); %average power across frequencies (so mean dB from ~40-80kHz, over temporal smoothing filter); do NOT normalize for now
% for filtering/smooothing:
wndsz = round(0.05/Tperiod); %convert seconds to samples
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
 usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 
AUCpower(j) = trapz(usvPowerPerSample)/totalSec; 
% AUCpowerARRAY(i) = AUCpower; 
%usvPwrPerSampleAll(end+1,1:22400) = usvPowerPerSample(1:22400); %ignore:
%this was written in case sampling is off between subjects
    
%change minpeakdistance depending on your expected duration of USVs 
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvs(j) = length(locUsvs); %number of USVs
numUsvs0to5(j)= length ( find (locUsvs<3750));
numUsvs5to10(j)= length ( find (locUsvs<7500 & locUsvs >=3750));
numUsvs10to15(j)= length ( find (locUsvs<11250 & locUsvs >=7500));
numUsvs15to20(j)= length ( find (locUsvs<15000 & locUsvs >=11250));
numUsvs20to25(j)= length ( find (locUsvs<18750 & locUsvs >=15000));
avgPeaks(j)= mean(pks); %average power of USVs per trial
stdPeaks(j) =std(pks); %standard deviation of USV power per trial

timepointUSVs= T(locUsvs); %timepoints that the USVs occurred at
if locUsvs >1 
%     timepoint= find (timepointUSVs>5,1);
%     firstcall= timepointUSVs(timepoint);
firstcall= timepointUSVs(1);
%     if latencytocall(j)< 0 % sometimes the noise is recognized as USV, in this case use the first call after 5s
%     firstcall= timepointUSVs(2); % only set to 2nd call, might need to change depend on data
%     end
else
    firstcall= NaN ; 
    latencytocall(j) = NaN;
end
latencytocall(j) = firstcall-5; % in seconds
interUSVinterval= diff(timepointUSVs); %find the difference between timepoints 
avgIUI(j) = mean(interUSVinterval);  %average the inter USV interval per trial
stdIUI(j) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial
    
if plotandsaveindividual     
disp(['Plotting...'])
% lowFreq = 1; 
% highFreq = length(F);
fig=figure;

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
imagesc(T, F(lowFreq:highFreq)./1000, signalCleaned(lowFreq:highFreq,:)); 
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
ylabel('USV power (total dB in 40-80kHz band)');
xlabel('time (seconds)')

%add by Jingyi, save fig

filname= [files(j).name];
saveas(fig, [wavPath, filname, 'USVPower.tif'], 'tif');
 
%now save individual data into matfile
matName = [wavPath, filname, '_usvData.mat']; %save data for each individual mouse
save(matName,'T','validSamples','usvPowerPerSampleSmooth', 'timepointUSVs', 'interUSVinterval', 'pks');
end

usvPowerPerSampleSmoothMat = cat(1, usvPowerPerSampleSmoothMat,...
         usvPowerPerSampleSmooth); % PETE- originally concat usvPowerPerSampleSmooth

    fprintf(1, '\n');               % put a line break between each sample
    % saveas (fig,'tif');
end
USVStats.AUCpower = AUCpower;
USVStats.NumUSV = numUsvs;
USVStats.AvgPeaks =avgPeaks;
USVStats.stdPeaks=stdPeaks;
USVStats.latencytoUSV =latencytocall;
USVStats.AvgIUI =avgIUI;
USVStats.stdIUI =stdIUI;
USVStats.numUsvs0to5= numUsvs0to5;
USVStats.numUsvs5to10= numUsvs5to10;
USVStats.numUsvs10to15= numUsvs10to15;
USVStats.numUsvs15to20= numUsvs15to20;
USVStats.numUsvs20to25= numUsvs20to25;
        matName = [wavPath, 'allstats.mat']; %save data for all mouse in the same folder
        save(matName, 'USVStats');
% AUCAvg= mean(AUCpower); 
% % 
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
% %     

%% heatmap - pete with smell timing, !!! IF THIS is giving an alpha data error, just copy paste this section into command window and it'll work

% smellTiming = xlsread('E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\Jason-WT SMUFF\2019_9_Sniff_vs_USV_timesheet\2019_9_Sniff_Timesheet.xlsx');
% smellMatrix = smellTiming(2:end, 2:end); % get rid of labels and unecessary vals
% smellMatrix = horzcat(zeros(size(smellMatrix)), smellMatrix); % add the zeros data before stimulus
% smellMatrix = repelem(smellMatrix,1, 750);
% sortSmells = [8,9,10,11,12,5,6,7,1,2,3,4]; % matlab reads files in different order than excel sheet
% smellMatrix = smellMatrix(sortSmells,:);
% 
% red = cat(3, ones(size(usvPowerPerSampleSmoothMat)), zeros(size(usvPowerPerSampleSmoothMat)), zeros(size(usvPowerPerSampleSmoothMat)));
% imagesc(usvPowerPerSampleSmoothMat);
% colormap(flipud(gray));
% hold on;
% h = imagesc(red, 'AlphaData', 1);
% hold off;
% set(h, 'AlphaData', smellMatrix(:,1:end-37)*0.4)
% colorbar;
% ax = gca;
% ax.XTick = (0:length(usvPowerPerSampleSmoothMat)/4:length(usvPowerPerSampleSmoothMat));
% ax.XTickLabel = ({'4', '1', '2', '3'});

% %% Plot the average signal -ZG
% 
% disp(['Loop finished. Now plotting average USV power...']);
% 
% % delete first row of zeros or the average suffers
% usvPowerPerSampleSmoothMat = usvPowerPerSampleSmoothMat(2:end,:); 
% 
% % take the average along the columns
% usvPowerPerSampleSmoothAvg = mean(usvPowerPerSampleSmoothMat);
% 
% %apply another round of filtering/smoothing
% wndsz = round(0.05/Tperiod);   % convert seconds to samples, original was 0.05
% gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
% gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
% validSamples =  length(usvPowerPerSample)-wndsz+1; % - PETE, 22499 if not smoothing
% usvPowerPerSampleSmoothAvg2 = conv(usvPowerPerSampleSmoothAvg, gaussFilter, 'same'); % not used - PETE
% 

% plot the average of the smoothed out signals
figure;
hold on;
% boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), std(usvPowerPerSampleSmoothMat(:, :, 1)),'-r');%red color for female
% boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), std(usvPowerPerSampleSmoothMat(:, :, 1))/sqrt(size(usvPowerPerSampleSmoothMat, 1)),'-m');%Magenda for GAD CNO days
% boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), 1.96*std(usvPowerPerSampleSmoothMat(:, :, 1))/sqrt(size(usvPowerPerSampleSmoothMat, 1)),'-b');%blue color for male or saline day for DREADD
boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), 1.96*std(usvPowerPerSampleSmoothMat(:, :, 1))/sqrt(size(usvPowerPerSampleSmoothMat, 1)),'-r');%red color for female or WT female cue
% orange= [ 0.9100 0.4100 0.1700];
% boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat), 1.96*std(usvPowerPerSampleSmoothMat)/sqrt(size(usvPowerPerSampleSmoothMat, 1)),'cmap', orange);%orange color for DREADD GAQ
% boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), std(usvPowerPerSampleSmoothAvg(:, :, 1)),'-k');%black color
% plot(T(1:validSamples), usvPowerPerSampleSmoothAvg, 'Color', [ 0.5843 0.8157 0.9882]); % r: - PETE
title(''); % note: next one overlaps it anyway
ylabel({'USV power/frequency';'(dB/Hz)'});
xlabel('time (seconds)');
%plot the stimulus box
xLims = get(gca, 'XLim');
axis([xLims(1) xLims(2) -0.15 4]);    % -0.05 to go slightly below lowest call, 3 should be above highest call for all VGlu T1T2 ChR2 and female ChR2
% axis([xLims(1) xLims(2) -0.15 1.2]);     % -0.05 to go slightly below lowest call, 2 should be above highest call for all other ChR2
% axis([xLims(1) xLims(2) -0.15 3.5]);  % -0.05 to go slightly below lowest call, 2.5 should be above highest call for all VGlu T3 ChR2
x = [5 10 10 5];
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(1) yLims(2) yLims(2)];
% patch(x, y, [1 0.3 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % orange box
% patch(x, y, [0 0.4 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % blue BOX FOR STIMULATION PERIOD
set(gca, 'FontSize', 20);     % font size of 20
hold off;
% 
% %PLOT heatmap with latency to call sorted --Add by JC on 12/2/2019
% Arowmax = max(usvPowerPerSampleSmoothMat, [], 2);
% [~,idx] = sort(Arowmax, 'descend');
% SortedMatrix = usvPowerPerSampleSmoothMat(idx,:);
% SortedNameList = nameList(idx); % update the index as well
% usvPowerPerSampleSmoothMat = SortedMatrix;
% 
% 
% %%report calling trials
% % length (USVStats.NumUSV);
% % length (find (USVStats.NumUSV>2))
% % length (find (USVStats.NumUSV>3))