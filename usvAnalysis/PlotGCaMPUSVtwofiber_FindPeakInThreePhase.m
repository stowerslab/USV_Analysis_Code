clear all; close all;
%%%Plot GCaMP with USV power, with behavior annotation
wavPath =  '/Users/jingyi/Dropbox (Scripps Research)/Jingyi Data backup/GCaMP Analysis/MasterDataset/';
folderPath =  '/Users/jingyi/Dropbox (Scripps Research)/Jingyi Data backup/GCaMP Analysis/MasterDataset/';

fontSz = 20;
mouseNums = {'BG6_T6R'};
totalMice = size(mouseNums, 2); 

gcampFrameRate= 20;
GCaMPtimesCell = {}; % photometry data timebase
frameTimesCell = {};

%extended time (s) from onset to draw figure
timeExtension = 1;
%a local time (s) window from onset to search for peaks
localWindow = 0.25;
%store processed data
avg_sig_norm_cleaned = zeros(totalMice,length(-timeExtension:1/gcampFrameRate:timeExtension));
avg_scrambledSig_norm = zeros(totalMice,length(-timeExtension:1/gcampFrameRate:timeExtension));
vAllOffset = 1;%offset into video parameter spreadsheet; for syncing to video only
vSheet = [wavPath 'MasterTimeGamp'];% read masterdataset
vData = xlsread(vSheet);

findPeaksFU=0; % find peaks in female urine window
findPeaksFemale =1; % find peaks in female window
findPeaksSpon =1; % find peaks in spontaneous window

for k = 1:totalMice 
    %read normalized GCamP data processed using Analyze_FIP_gui
    mouseNum = mouseNums{1,k};
    matName = [mouseNum, '_processed.mat'];
    matName = char(matName);
    load(fullfile(folderPath,matName)); 
    sig_norm = sig_norm * 100; % calculate %
    GCaMPCell{k}=sig_norm;
    numGcampSamples = length(GCaMPCell{k});
    totalGCampSec= length(sig_norm)/gcampFrameRate;
    timescaleGcamp= 0:1/gcampFrameRate:totalGCampSec-(1/gcampFrameRate);
%     GCaMPtimesCell{k} = 0:gcampFrameRate:(numGcampSamples*gcampFrameRate)-gcampFrameRate; % 
    
AUCpowerARRAY = [];
usvPwrPerSampleAll = [];    
%Read USV file that has been trimed to match GCamP recording time       
waveFile = [mouseNum,'_FP.wav'];
disp(['Reading WAV file: ', waveFile])
 [y, Fs] = audioread(fullfile(folderPath,waveFile));
i = audioinfo(fullfile(folderPath,waveFile));
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

refPower = 10^-12; %reference power is 10?12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
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

%find syllables with minimum distance
minSyllableDistance = 2; %define a sentence by 2s breaks
syllable = zeros(length(timescaleUSV),1);
syllable(locUsvs,1)=1;
[~,locSyllable] = findpeaks(syllable,'MinPeakDistance',minSyllableDistance/Tperiod);

avg_sig_norm = zeros(length(locSyllable),length(-timeExtension:1/gcampFrameRate:timeExtension));

for n = 1:length(locSyllable)
    %preid is the closest time point in gcamp for each onset
    [~,preid] = min(abs(timescaleGcamp-locSyllable(n)*Tperiod));
    %find the most prominent peak within a local time window to correct misalignment
    [~,locOnset,~,prm] = findpeaks(sig_norm(1,preid-localWindow*gcampFrameRate:preid+localWindow*gcampFrameRate));
    if length(locOnset)>=1
        %correct misalignment
        postid = preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 ;
        avg_sig_norm(n,:) =sig_norm(1,postid-timeExtension*gcampFrameRate:postid+timeExtension*gcampFrameRate);
    end
end
avg_sig_norm_cleaned(k,:) = mean(avg_sig_norm(avg_sig_norm(:,1)~=0,:));

numScrambled = size(avg_sig_norm(avg_sig_norm(:,1)~=0,:),1); % choose the same number of scrambled sample as USV sentence onset detected in the file
scrambledID =  datasample(1:length(sig_norm(1,1:numGcampSamples-timeExtension*gcampFrameRate)),numScrambled);
scrambledSig_norm = zeros(numScrambled,length(-timeExtension:1/gcampFrameRate:timeExtension));
for j = 1:numScrambled
    scrambledSig_norm(j,:) = sig_norm(1,scrambledID(j)-timeExtension*gcampFrameRate:scrambledID(j)+timeExtension*gcampFrameRate);
end
avg_scrambledSig_norm(k,:) = mean(scrambledSig_norm);

figure;% plot for the entire data set
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
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([0,totalGCampSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
subplot(3,1,3);
plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned(k,:),'green');
hold on
plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_scrambledSig_norm(k,:),'black');
title('GCaMP signals and Scrambled signal at onset');
ylabel ('dF/F (%)');
xlabel ('time (s)')
% subplot(4,1,4);
% plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_scrambledSig_norm(k,:),'black');
% title('scrambled GCaMP signals');
% ylabel ('dF/F (%)');
% xlabel ('time (s)')


if findPeaksFU
 syncVidSecond = vData(vAllOffset+k,1);
 stimStartVidSec = vData(vAllOffset+k,2) - syncVidSecond; %subtracted sync time from start time
 stimEndVidSec = vData(vAllOffset+k,3) - syncVidSecond-5; % subtracted sync time from end time and substract 5s for opening the cage for next stimulation
 locSyllableFUId = find (locSyllable >stimStartVidSec/Tperiod & locSyllable < stimEndVidSec/Tperiod);
 locSyllableFU = locSyllable(locSyllableFUId);
 avg_sig_norm_FU = zeros(length(locSyllableFU),length(-timeExtension:1/gcampFrameRate:timeExtension));

for m = 1:length(locSyllableFU)
    %preid is the closest time point in gcamp for each onset
    [~,preid] = min(abs(timescaleGcamp-locSyllableFU(m)*Tperiod));
    %find the most prominent peak within a local time window to correct misalignment
    [~,locOnset,~,prm] = findpeaks(sig_norm(1,preid-localWindow*gcampFrameRate:preid+localWindow*gcampFrameRate));
    if length(locOnset)>=1
        %correct misalignment
        postid = preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 ;
        avg_sig_norm_FU(m,:) =sig_norm(1,postid-timeExtension*gcampFrameRate:postid+timeExtension*gcampFrameRate);
    end
end
avg_sig_norm_cleaned_FU(k,:) = mean(avg_sig_norm_FU(avg_sig_norm_FU(:,1)~=0,:));


% Now plot data 
figure;% plot for the female part only
subplot(3,1,1);
plot(timescaleUSV, usvPowerPerSampleSmooth, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([stimStartVidSec,stimEndVidSec]);
title('USV power');
% xlim([0,totalGCampSec]);
ylabel ('AUC power (dB)');
xlabel ('time (s)')
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([stimStartVidSec,stimEndVidSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)');
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([stimStartVidSec,stimEndVidSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
subplot(3,1,3);
plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned_FU(k,:),'green');
title('GCaMP signals at onset during female Urine');
ylabel ('dF/F (%)');
xlabel ('time (s)')

end


if findPeaksFemale
 syncVidSecond = vData(vAllOffset+k,1);
 stimStartVidSec = vData(vAllOffset+k,3) - syncVidSecond; %subtracted sync time from start time
 stimEndVidSec = vData(vAllOffset+k,4) - syncVidSecond-5; % subtracted sync time from end time and substract 5s for opening the cage for next stimulation
 
 locSyllableFId = find (locSyllable >stimStartVidSec/Tperiod & locSyllable < stimEndVidSec/Tperiod);
 locSyllableF = locSyllable(locSyllableFId);
 avg_sig_norm_F = zeros(length(locSyllableF),length(-timeExtension:1/gcampFrameRate:timeExtension));

for m = 1:length(locSyllableF)
    %preid is the closest time point in gcamp for each onset
    [~,preid] = min(abs(timescaleGcamp-locSyllableF(m)*Tperiod));
    %find the most prominent peak within a local time window to correct misalignment
    [~,locOnset,~,prm] = findpeaks(sig_norm(1,preid-localWindow*gcampFrameRate:preid+localWindow*gcampFrameRate));
    if length(locOnset)>=1
        %correct misalignment
        postid = preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 ;
        avg_sig_norm_F(m,:) =sig_norm(1,postid-timeExtension*gcampFrameRate:postid+timeExtension*gcampFrameRate);
    end
end
avg_sig_norm_cleaned_F(k,:) = mean(avg_sig_norm_F(avg_sig_norm_F(:,1)~=0,:));


% Now plot data 
figure;% plot for the female part only
subplot(3,1,1);
plot(timescaleUSV, usvPowerPerSampleSmooth, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([stimStartVidSec,stimEndVidSec]);
title('USV power');
% xlim([0,totalGCampSec]);
ylabel ('AUC power (dB)');
xlabel ('time (s)')
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([stimStartVidSec,stimEndVidSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)');
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([stimStartVidSec,stimEndVidSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
subplot(3,1,3);
plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned_F(k,:),'green');
title('GCaMP signals at onset during female period');
ylabel ('dF/F (%)');
xlabel ('time (s)')
end

if findPeaksSpon
 syncVidSecond = vData(vAllOffset+k,1);
 stimStartVidSec = vData(vAllOffset+k,4) - syncVidSecond; %subtracted sync time from start time
 stimEndVidSec = vData(vAllOffset+k,5) - syncVidSecond; % subtracted sync time from end time 
 
locSyllableSId = find (locSyllable >stimStartVidSec/Tperiod & locSyllable < stimEndVidSec/Tperiod);
 locSyllableS = locSyllable(locSyllableSId);
 avg_sig_norm_S = zeros(length(locSyllableS),length(-timeExtension:1/gcampFrameRate:timeExtension));

for m = 1:length(locSyllableS)
    %preid is the closest time point in gcamp for each onset
    [~,preid] = min(abs(timescaleGcamp-locSyllableS(m)*Tperiod));
    %find the most prominent peak within a local time window to correct misalignment
    [~,locOnset,~,prm] = findpeaks(sig_norm(1,preid-localWindow*gcampFrameRate:preid+localWindow*gcampFrameRate));
    if length(locOnset)>=1
        %correct misalignment
        postid = preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 ;
        avg_sig_norm_S(m,:) =sig_norm(1,postid-timeExtension*gcampFrameRate:postid+timeExtension*gcampFrameRate);
    end
end
avg_sig_norm_cleaned_S(k,:) = mean(avg_sig_norm_S(avg_sig_norm_S(:,1)~=0,:));


% Now plot data 
figure;% plot for the female part only
subplot(3,1,1);
plot(timescaleUSV, usvPowerPerSampleSmooth, 'r');% plot USV power
hold on;
plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'b*')
xlim([stimStartVidSec,stimEndVidSec]);
title('USV power');
% xlim([0,totalGCampSec]);
ylabel ('AUC power (dB)');
xlabel ('time (s)')
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
xlim([stimStartVidSec,stimEndVidSec]);
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)');
subplot(3,1,2);
plot(timescaleGcamp, sig_norm(1,:), 'g'); % plot GCaMP Trace
hold on;
plot(locSyllable*Tperiod, ones(length(locSyllable),1).*max(sig_norm(1,:)), 'b+')
xlim([stimStartVidSec,stimEndVidSec]);
title('GCaMP signals');
ylabel ('dF/F (%)');
xlabel ('time (s)')
subplot(3,1,3);
plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned_S(k,:),'green');
title('GCaMP signals at onset during Spontaneous calling');
ylabel ('dF/F (%)');
xlabel ('time (s)')
end
end

