% top-level script to process USV audio data from pup isolation plot power is USV band
% NOTE: requires Statistics toolbox

clear all; close all;
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_11_pupUsvs\cd1_p7\';
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_12_pups\b6p14\';
wavPath = 'F:\StowersLabData\Jason\data\behaviorVideoUsv\2018_12_pupUsvs\b6p14\';
% wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_09_pupUsvs\p7_pretest\';

% mouseNums = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'};
% mouseNums = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14','15','16'};
% mouseNums = {'4', '5', '9'};
mouseNums = {'SeaWave _20181205 151749'};
totalMice = size(mouseNums, 2);

numUsvs = zeros(totalMice,1);

for i = 1:totalMice
    mouseNum = mouseNums{1,i};
    waveFile = [wavPath, mouseNum, '.wav'];
    disp(['Reading WAV file: ', waveFile])
    [y, Fs] = audioread(waveFile);
    Yinfo= audioinfo(waveFile);
    totalSec =Yinfo.Duration;
    TwoMinSamples = Fs*120;

    % trim file if necessary by % of total length or seconds
%     startSec = 120;
%     stopSec =125;
%     startPercent = startSec/totalSec;
%     stopPercent = stopSec/totalSec;
%     startInd = round(length(y)*startPercent);
%     stopInd = round(length(y)*stopPercent);
    startInd = 1;
    stopInd = length(y);
%     startInd = stopInd - TwoMinSamples;
    yTrim = y(startInd:stopInd); 
    calcDuration = (stopInd-startInd)/Fs;

    disp(['Taking FFT...'])
    nfft = 512;
    window = 512;
    noverlap = window*0.5;
    % thresh = -90; % threshold in decibels;
    [~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
    %note that P is the power spectral density in W/Hz
    Tperiod =calcDuration/length(T);

    refPower = 10^-12; %eference power is 10^-12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
    signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 

    disp(['Filtering noise...'])
    % idea here take z-score & compare to other freqs to remove broadband noise: 
    zsignal = zscore(signal);
    lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
    highFreq = find(F>100000,1,'first'); % index for high cutoff used below
    zsignal(1:lowFreq,:) = 0; % lowpass - set everything below cutoff to zero
    zthresh = 1.3; %rahter arbitrary
    zsignal(zsignal<zthresh) = 0; %threshold zscore
    signalCleaned = signal; %create a copy to clean below
    signalCleaned(zsignal==0) = 0; %JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!

    % calculate power in the whistle snippets over time
    disp(['Calculating acoustic power...'])
    usvPowerPerSample = mean(signalCleaned); %average power across frequencies (so mean dB/Hz from lowFreq-highFreq, over temporal smoothing filter); do NOT normalize for now
    % for filtering/smooothing:
    wndSec = 0.1;
    wndsz = round(wndSec/Tperiod); %convert seconds to samples
    gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
    gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
    validSamples = length(usvPowerPerSample)-wndsz+1;
     usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 

     %count USVs based on power peaks above threshold
     powerThresh = 1; 
    [pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
    numUsvs(i) = length(locUsvs);
    
    figure;
    subplot(2,1,1)
    title('original');
    imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
    set(gca,'YDir','normal') %flip to make small freq on bottom
    colormap(1-bone);
    ylabel('frequency (kHz)');
    xlabel('time (seconds)')

    % subplot(3,1,2)
    % title('filtered, z-scored, and thresholded');
    % imagesc(T, F(lowFreq:highFreq)./1000, signalCleaned(lowFreq:highFreq,:)); 
    % set(gca,'YDir','normal') %flip to make small freq on bottom
    % colormap(1-bone);
    % ylabel('frequency (kHz)');
    % xlabel('time (seconds)')

    subplot(2,1,2)
    title('power per sample');
    plot(T(1:validSamples), usvPowerPerSampleSmooth)
    hold on;
    plot([0,calcDuration], [powerThresh, powerThresh], 'r--')
    plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'g*')
    ylabel('USV power (mean dB/Hz in 40-100kHz band)');
    xlabel('time (seconds)')
    hold off;

end

figure;
plot(numUsvs, 'k*')
ylabel('# USVs');
xlabel('mouse #')







