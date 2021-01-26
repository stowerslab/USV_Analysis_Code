%JC for sentence analysis
%break USV files to sentences and calculate Number of USVs/sentence;
%sentence length etc for group comparisons
clear all; close all;
wavPath = '\\172.29.164.29\jingyi chen\Dropbox\Jingyi Data backup\Behavior data\WT behaviors\2018-2 WT USVs\2018-2-26 wt USV FEMALE\2MINS TRIMED AUDIO\';
%Note have to put allstats.mat into another foder
files = dir('*.wav');
totalMice = size(files, 1);
usvPwrPerSampleAll = [];
timepointUSVs = [];
nameList = {};
usvPowerPerSampleSmoothMat = zeros(1,89963); % for 2 mins WT SMUF analysis only

for j = 1:totalMice

waveFile = [wavPath, files(j).name];
nameList = [nameList files(j).name];
disp(['Reading WAV file: ', waveFile])
[y, Fs] = audioread(waveFile);

% % %  error checking for different sample rates -ZG
% %     disp(['Checking the sampling rate: ', num2str(Fs)])
% %     if Fs == 250000
% %         y = resample(y, 192000, 250000);
% %     end
% %     Fs = 192000;
    
i = audioinfo(waveFile);
totalSec = i.Duration;

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
usvPowerPerSample = mean(signalCleaned); %average power across frequencies (so mean dB from ~40-80kHz, over temporal smoothing filter); do NOT normalize for now
% for filtering/smooothing:
wndsz = round(0.05/Tperiod); %convert seconds to samples
gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
validSamples = length(usvPowerPerSample)-wndsz+1;
 usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 
% AUCpower(j) = trapz(usvPowerPerSample)/totalSec; 
% AUCpowerARRAY(i) = AUCpower; 
%usvPwrPerSampleAll(end+1,1:22400) = usvPowerPerSample(1:22400); %ignore:
%this was written in case sampling is off between subjects
    
%change minpeakdistance depending on your expected duration of USVs 
powerThresh = 1; 
[pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
numUsvs(j) = length(locUsvs); %number of USVs
avgPeaks(j)= mean(pks); %average power of USVs per trial
stdPeaks(j) =std(pks); %standard deviation of USV power per trial

timepointUSVs= T(locUsvs); %timepoints that the USVs occurred at

interUSVinterval= diff(timepointUSVs); %find the difference between timepoints 
avgIUI(j) = mean(interUSVinterval);  %average the inter USV interval per trial
stdIUI(j) = std(interUSVinterval);  %standard deviation of the interUSV interval per trial


idx_end_sentences = []; % record the end of sentences
idx_init_sentences = []; % record the begining of sentences
manipulate_timepoint = timepointUSVs;
for i=1:length(interUSVinterval)
%     if interUSVinterval(i)>=2 && interUSVinterval(i+1) >= 2 % set the gap between sentences is 2s
%         manipulate_timepoint(i) = 0;
%         manipulate_timepoint(i+1) = 0;
%     elseif interUSVinterval(i)>=2
    if interUSVinterval(i)>=2
        idx_end_sentences = [idx_end_sentences,i];
        idx_init_sentences = [idx_init_sentences,i+1];
    end
end
idx_end_sentences = [idx_end_sentences,length(timepointUSVs)];

for i=1:length(idx_init_sentences)
    end_sen = idx_end_sentences(i+1); % extract the idx of last peaks in a sentence
    init_sen = idx_init_sentences(i);
%     disp(['init_sen: ', init_sen]);
    gap = timepointUSVs(end_sen)-timepointUSVs(init_sen);
%     disp(gap);
    if gap < 1
%         disp(timepointUSVs(end_sen)-timepointUSVs(init_sen));
%         disp('in the loop');
      idx_end_sentences(i+1)=0;
      idx_init_sentences(i)=0;
     end  
%     disp(['ith round: ',i]);
%     disp(['init: ', idx_init_sentences(i)]);
%     disp(['end: ', idx_end_sentences(i+1)]);
end

idx_end_sentences = idx_end_sentences(idx_end_sentences ~= 0);
idx_init_sentences = idx_init_sentences(idx_init_sentences ~= 0);

%Now extract sentence feature
for k=1:length(idx_init_sentences)
    end_sen = idx_end_sentences(k+1); % extract the idx of last peaks in a sentence
    %disp(['end_sen: ', end_sen]);
    init_sen = idx_init_sentences(k);
    %disp(['init_sen: ', init_sen]);
    onset = locUsvs(init_sen)+20;
    offset = locUsvs(end_sen)+20;
    sen_usv = usvPowerPerSampleSmooth(onset:offset);
    %save(sen_file_usv,'sen_usv','-append');
    [Sen_pks, Sen_locUsvs] = findpeaks(sen_usv,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %run finding peaks again with trimed sentence
    Sen_numUsvs(k) = length(Sen_locUsvs); %number of USVs
    Sen_avgPeaks(k)= mean(Sen_pks); %average power of sentence USVs per trial
    Sen_stdPeaks(k) =std(Sen_pks); %standard deviation of sentence USV power per trial
    init_sen = idx_init_sentences(k); % update the idx of the next peaks
   
end
filname= [files(j).name];
matName = [wavPath, filname, '_SentenceData.mat']; %save data for each individual mouse
save(matName,'Sen_numUsvs','Sen_avgPeaks', 'Sen_stdPeaks');
end



