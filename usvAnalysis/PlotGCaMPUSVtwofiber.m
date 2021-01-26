clear all; close all;
%%%Plot GCaMP with USV power, with behavior annotation
wavPath =  '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\GCaMP Analysis\2018-12 GCaMP Batch 3\BG9\';

fontSz = 20;
mouseNums = {'BG14_T2'};
totalMice = size(mouseNums, 2); 
plotSeparate = 1;

gcampFrameRate= 20;
GCaMPtimesCell = {}; % photometry data timebase
frameTimesCell = {};
vmousNums = char(mouseNums);

%%plot each annotated behavior with GCaMP


for k = 1:totalMice %#ok<UNRCH>
    %read normalized GCamP data processed using Analyze_FIP_gui
    mouseNum = mouseNums{1,k};
    matName = [mouseNum, '_processed.mat'];
    matName = char(matName);
    load(matName); 
%     sig_norm = sig_norm * 100 % converted to %
    sig_norm = sig_norm * 100; % calculate %
    GCaMPCell{k}=sig_norm;
    numGcampSamples = length(GCaMPCell{k});
    totalGCampSec= length(sig_norm)/gcampFrameRate;
    timescaleGcamp= 0:1/gcampFrameRate:totalGCampSec-(1/gcampFrameRate);
%     GCaMPtimesCell{k} = 0:gcampFrameRate:(numGcampSamples*gcampFrameRate)-gcampFrameRate; % 
    
AUCpowerARRAY = [];
usvPwrPerSampleAll = [];    
%Read USV file that has been trimed to match GCamP recording time       
waveFile = [vmousNums,'_FP.wav'];
disp(['Reading WAV file: ', waveFile])
 [y, Fs] = audioread(waveFile);
i = audioinfo(waveFile);
totalSec = i.Duration;
USVSamplerate=i.SampleRate;

% becuase the USV file has already been trimed so use following code
startInd = 1;
stopInd = length(sig_norm)/gcampFrameRate * USVSamplerate; % to match time, should be less than 0.5 s different
yTrim = y(startInd:stopInd);

disp(['Taking FFT...'])
nfft = 512;
window = 512;
noverlap = window*0.5;
% thresh = -90; % threshold in decibels;
[~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
%note that P is the power spectral density in W/Hz
Tperiod =i.Duration/length(T);

refPower = 10^-12; %eference power is 10?12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 

disp(['Filtering noise...'])
% idea here take z-score & compare to other freqs to remove broadband noise: 
zsignal = zscore(signal);
lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
highFreq = find(F>80000,1,'first'); % index for high cutoff used below
zsignal(1:lowFreq,:) = 0; % lowpass - set everything below cutoff to zero
zsignal(highFreq:end,:) = 0; % highpass - add by Jingyi on 11/29/2018
zthresh = 1.5; %1.3 for pup USV
zsignal(zsignal<zthresh) = 0; %threshold zscore
signalCleaned = signal; %create a copy to clean below
signalCleaned(zsignal==0) = 0; %JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

% calculate power in the whistle snippets over time
disp(['Calculating acoustic power...'])
usvPowerPerSample = mean(signalCleaned); %average power across frequencies (so mean dB from ~40-80kHz, over temporal smoothing filter); do NOT normalize for now
% for filtering/smooothing:
wndsz = round(0.05/Tperiod); %convert seconds to samples
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
 usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 
AUCpower = trapz(usvPowerPerSample)/totalSec; 

% usvPwrPerSampleAll(end+1,1:22400) = usvPowerPerSample(1:22400);
samplerateRatio= length (usvPowerPerSampleSmooth)/length(sig_norm);
timescaleUSV=0:1/(gcampFrameRate*samplerateRatio):totalGCampSec-1/(gcampFrameRate*samplerateRatio); % match USV time with GCamP time

AUCpower(k) = trapz(usvPowerPerSample)/totalSec; 
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvs(k) = length(locUsvs); %number of USVs

% avgPeaks(k)= mean(pks); %average power of USVs per trial
% stdPeaks(k) =std(pks); %standard deviation of USV power per trial

% timepointUSVs= T(locUsvs); %timepoints that the USVs occurred at
% if locUsvs >1 
%     firstcall= timepointUSVs(1); 
% else
%     firstcall= NaN ; 
% end
% latencytocall(k) = firstcall-5; %(in seconds) 
% interUSVinterval= diff(timepointUSVs); %find the difference between timepoints 
% avgIUI(k) = mean(interUSVinterval);  %average the inter USV interval per trial
% stdIUI(k) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial


if plotSeparate
% sig_norm_down = downsample(sig_norm,20);%downsample GCamP to match behavior annotation happened in seconds
figure;
subplot(3,1,1);
plot(timescaleUSV, usvPowerPerSampleSmooth, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([0,totalGCampSec]);
title('USV power');
% xlim([0,totalGCampSec]);
ylabel ('AUC power (dB)');
xlabel ('time (s)')
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plor GCaMP Trace
xlim([0,totalGCampSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
ylim([-4 12]);
subplot(3,1,3);
plot(timescaleGcamp, sig_norm(2,:), 'g'); % plor GCaMP Trace
xlim([0,totalGCampSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
ylim([-4 12]);
% xlabel ('time (s)')
% subplot(3,1,3);% behavior annotation
% heatmap (annotation,'Colormap',parula(40),'MissingDataColor',[1 1 1],'GridVisible','off'); %use parula 40 to seperate color, missing data set to white
% title('behavior annotations');
% filname= [files(j).name];
% saveas(figure, [wavPath, mouseNums, '_GCaMP.tif'], 'tif');%save the figure


else
% %Plot GCamP together with USV power individually
cmapGcamp = [0 0.7 0];
cmapUSV = [1 0 0];
hold on;
xlabel('time (sec)', 'FontSize', fontSz);
ylabel('USV Power , dF/F', 'FontSize', fontSz);
plot(timescaleUSV, usvPowerPerSampleSmooth,  'Color', cmapUSV, 'LineWidth', 1)
ylim([-5 15]);
plot(timescaleGcamp, sig_norm(1,:), 'Color', cmapGcamp, 'LineWidth', 1);
ylim([-5 15]);
xlim([0,totalGCampSec]);
% axis tight;
hold off;


cmapGcamp = [0 0.7 0];
cmapUSV = [1 0 0];
hold on;
xlabel('time (sec)', 'FontSize', fontSz);
ylabel('USV Power , dF/F', 'FontSize', fontSz);
plot(timescaleUSV, usvPowerPerSampleSmooth,  'Color', cmapUSV, 'LineWidth', 1)
ylim([-5 15]);
plot(timescaleGcamp, sig_norm(2,:), 'Color', cmapGcamp, 'LineWidth', 1);
ylim([-5 15]);
xlim([0,totalGCampSec]);
% axis tight;
hold off;

end
end

