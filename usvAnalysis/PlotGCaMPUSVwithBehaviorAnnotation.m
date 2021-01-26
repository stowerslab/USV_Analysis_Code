clear all; close all;
%%%Plot GCaMP with USV power, with behavior annotation
wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2018-9 GCaMP BNST\BG6\';

fontSz = 20;
mouseNums = {'BG6_T4R'};
totalMice = size(mouseNums, 2); 
Femaleurine=1;
Female=1;
gcampFrameRate= 20;
% GCaMPtimesCell = {}; % photometry data timebase
% frameTimesCell = {};
% frameTimesCell{k} = frameTimes.*60; %convert back to seconds

%%plot each annotated behavior with GCaMP
vmousNums = char(mouseNums); 
if Femaleurine
vBehaviorFU = xlsread(vmousNums,'FU'); % read behavior annotation for FemaleUrine
annotationFU = vBehaviorFU (:,2);
annotationFU = annotationFU.' ;% transpose data for heatmap
vSyn=vBehaviorFU (1,4); % the second GCaMP starts
vStartFU=vBehaviorFU (1,3)- vSyn;% start second
vTimeFU= length (vBehaviorFU(:,3));% how many seconds for annotation
vEndFU =vStartFU + vTimeFU; % calculate ending seconds
 for k = 1:totalMice %#ok<UNRCH>
    %read normalized GCamP data processed using Analyze_FIP_gui
    mouseNum = mouseNums{1,k};
    matName = [mouseNum, '_processed.mat'];
    matName = char(matName);
    load(matName); 
    sig_norm = sig_norm * 100; % converted to %
    %%Cut GCaMP for Female Urine period
    GCaMPFU=sig_norm(1,vStartFU*gcampFrameRate:vEndFU*gcampFrameRate-1); 
    numGcampSamplesFU = length(GCaMPFU);
    totalGCampSecFU= numGcampSamplesFU/gcampFrameRate; % convert it to seconds just to check with vTime FU to make sure code works
    timescaleGcampFU= 0:1/gcampFrameRate:vTimeFU-1/gcampFrameRate; % convert to second
%     GCaMPtimesCell{k} = 0:gcampFrameRate:(numGcampSamples*gcampFrameRate)-gcampFrameRate; % 

AUCpowerFUARRAY = [];   
%Read USV file that has been trimed to match GCamP recording time       
waveFile = [vmousNums,'_FP.wav'];
disp(['Reading WAV file: ', waveFile])
 [y, Fs] = audioread(waveFile);
i = audioinfo(waveFile);
totalSec = i.Duration;
USVSamplerate=i.SampleRate;

% becuase the USV file has already been trimed so use following code to
% match GCaMP time scale, not precise timing but difference should be less
% than 0.5s

% Trim the data for FU and F seperately 
startIndFU = vStartFU * USVSamplerate;
stopIndFU = vEndFU * USVSamplerate;
yTrim = y(startIndFU:stopIndFU);

disp(['Taking FFT for Female Urine...'])
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
AUCpowerFU = trapz(usvPowerPerSample)/totalSec; 
AUCpowerFU(k) = trapz(usvPowerPerSample)/totalSec; 
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvsFU(k) = length(locUsvs); %number of USVs

avgPeaksFU(k)= mean(pks); %average power of USVs per trial
stdPeaksFU(k) =std(pks); %standard deviation of USV power per trial

timepointUSVsFU= T(locUsvs); %timepoints that the USVs occurred at
if locUsvs >1 
    firstcall= timepointUSVsFU(1); 
else
    firstcall= NaN ; 
end
latencytocallFU(k) = firstcall-5; %(in seconds) 
% interUSVintervalFU= diff(timepointUSVs); %find the difference between timepoints 
% avgIUIFU(k) = mean(interUSVinterval);  %average the inter USV interval per trial
% stdIUIFU(k) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial
usvPowerPerSampleFU = usvPowerPerSampleSmooth;
%Now trim the FU and F USV data using the smoothed data
% samplerateSmoothed = length (usvPowerPerSampleSmooth)/totalSec;
% 
% usvPowerPerSampleFU = usvPowerPerSampleSmooth (1, round(vStartFU*samplerateSmoothed):round(vEndFU*samplerateSmoothed));% trim the FU responding signals
% usvPowerPerSampleF = usvPowerPerSampleSmooth (1, round(vStartF*samplerateSmoothed):round(vEndF*samplerateSmoothed));% trim the F responding signals

samplerateRatioFU= length (usvPowerPerSampleFU)/length(GCaMPFU);
timescaleUSVFU=0:1/(gcampFrameRate*samplerateRatioFU):vTimeFU-1/(gcampFrameRate*samplerateRatioFU);  % match USV time with GCamP time


figure;
subplot(3,1,1);% behavior annotation
heatmap (annotationFU,'Colormap',spring(20),'MissingDataColor',[1 1 1],'GridVisible','off'); %use parula 40 to seperate color, missing data set to white
title('behavior annotations for female urine');
subplot(3,1,2);
plot(timescaleGcampFU, GCaMPFU, 'g'); % plor GCaMP Trace
xlim([0,vTimeFU]);
title('GCaMP signals for Female Urine');
ylabel ('dF/F (%)');
ylim([-4 12])
xlabel ('time (s)')
subplot(3,1,3);
plot(timescaleUSVFU, usvPowerPerSampleFU, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([0,vTimeFU]);
title('USV power');
ylabel ('AUC power (dB)');
xlabel ('time (s)')
ylim([0 5]);
%  saveas(figure, [wavPath, mouseNums, '_AnnotationFU.tif'], 'tif');%save the figure
 end
end

if Female
vBehaviorF = xlsread(vmousNums,'F'); % read behavior annotation for Female
% vBehavior(isnan(vBehavior))=0; % Set all NAN to 0
annotationF = vBehaviorF (:,2);
annotationF = annotationF.' ;% transpose data for heatmap
vSyn=vBehaviorF (1,4);
vStartF = vBehaviorF (1,3)-vSyn;% start second
vTimeF= length (vBehaviorF(:,3));% how many seconds for annotation
vEndF= vStartF+ vTimeF; 
for j = 1:totalMice %#ok<UNRCH>
    %read normalized GCamP data processed using Analyze_FIP_gui
    mouseNum = mouseNums{1,j};
    matName = [mouseNum, '_processed.mat'];
    matName = char(matName);
    load(matName); 
    sig_norm = sig_norm * 100; % converted to %
    %%Cut GCaMP for Female Urine period
    GCaMPF=sig_norm(1,vStartF*gcampFrameRate:vEndF*gcampFrameRate-1); 
    numGcampSamplesF = length(GCaMPF);
    totalGCampSecF= numGcampSamplesF/gcampFrameRate; % convert it to seconds just to check with vTime FU to make sure code works
    timescaleGcampF= 0:1/gcampFrameRate:vTimeF-1/gcampFrameRate; % convert to second
%     GCaMPtimesCell{k} = 0:gcampFrameRate:(numGcampSamples*gcampFrameRate)-gcampFrameRate; % 

AUCpowerFARRAY = [];   
%Read USV file that has been trimed to match GCamP recording time       
waveFile = [vmousNums,'_FP.wav'];
disp(['Reading WAV file: ', waveFile])
 [y, Fs] = audioread(waveFile);
i = audioinfo(waveFile);
totalSec = i.Duration;
USVSamplerate=i.SampleRate;

% becuase the USV file has already been trimed so use following code to
% match GCaMP time scale, not precise timing but difference should be less
% than 0.5s

% Trim the data for FU and F seperately 
startIndF = vStartF * USVSamplerate;
stopIndF = vEndF * USVSamplerate;
yTrim = y(startIndF:stopIndF);

disp(['Taking FFT for Female stimuli...'])
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
AUCpowerF = trapz(usvPowerPerSample)/totalSec; 
AUCpowerF(j) = trapz(usvPowerPerSample)/totalSec; 
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvsFU(j) = length(locUsvs); %number of USVs

avgPeaksF(j)= mean(pks); %average power of USVs per trial
stdPeaksF(j) =std(pks); %standard deviation of USV power per trial

timepointUSVsF= T(locUsvs); %timepoints that the USVs occurred at
if locUsvs >1 
    firstcall= timepointUSVsF(1); 
else
    firstcall= NaN ; 
end
latencytocallF(j) = firstcall-5; %(in seconds) 
% interUSVintervalFU= diff(timepointUSVs); %find the difference between timepoints 
% avgIUIFU(k) = mean(interUSVinterval);  %average the inter USV interval per trial
% stdIUIFU(k) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial
usvPowerPerSampleF = usvPowerPerSampleSmooth;
%Now trim the FU and F USV data using the smoothed data
% samplerateSmoothed = length (usvPowerPerSampleSmooth)/totalSec;
% 
% usvPowerPerSampleFU = usvPowerPerSampleSmooth (1, round(vStartFU*samplerateSmoothed):round(vEndFU*samplerateSmoothed));% trim the FU responding signals
% usvPowerPerSampleF = usvPowerPerSampleSmooth (1, round(vStartF*samplerateSmoothed):round(vEndF*samplerateSmoothed));% trim the F responding signals

samplerateRatioF= length (usvPowerPerSampleF)/length(GCaMPF);
timescaleUSVF=0:1/(gcampFrameRate*samplerateRatioF):vTimeF-1/(gcampFrameRate*samplerateRatioF);  % match USV time with GCamP time


figure;
subplot(3,1,1);% behavior annotation
heatmap (annotationF,'Colormap',parula(20),'MissingDataColor',[1 1 1],'GridVisible','off'); %use parula 40 to seperate color, missing data set to white
title('behavior annotations for female urine');
subplot(3,1,2);
plot(timescaleGcampF, GCaMPF, 'g'); % plor GCaMP Trace
xlim([0,vTimeF]);
title('GCaMP signals for Female Urine');
ylabel ('dF/F (%)');
ylim([-4 12])
xlabel ('time (s)')
subplot(3,1,3);
plot(timescaleUSVF, usvPowerPerSampleF, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([0,vTimeF]);
title('USV power');
ylabel ('AUC power (dB)');
xlabel ('time (s)')
ylim([0 5]);
%  saveas(figure, [wavPath, mouseNums, '_AnnotationF.tif'], 'tif');%save the figure
end
end
