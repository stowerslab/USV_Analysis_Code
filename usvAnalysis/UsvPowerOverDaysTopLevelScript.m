% top-level script to process USV audio data from pup isolation plot power is USV band
% NOTE: requires Statistics toolbox

clear all; close all;
wavPath = 'C:\Users\jason\OneDrive\stowersLab\data\behaviorVideoUsv\2018_12_pups\cd1\';

totalMice = 19;
%create filenames:
for i=1:totalMice
    prefix = 'p12_';
    mouseNums{1,i} = [prefix num2str(i)]; 
end

male = {'1','5','8','9','10','13','14','15','16','18'}; % else female

% declare variables to compute and save
numUsvs = zeros(totalMice,1);
usvPwrSum = zeros(totalMice,1);

% program control:
findUsvs =0; %first: manual step; do this only
combineMice = 0; %make groups
computeStats = 1; %compute and plot
normalize = 1; %normalize to max. USVs over multiple days?

matNameSaveSummary = [wavPath, 'allCD1_p12.mat'];

groupsToCompare = {'allCD1_p7', 'allCD1_p8', 'allCD1_p9', 'allCD1_p10', 'allCD1_p11', 'allCD1_p12'}; 

stat1 = 'numUsvs';
usvStats = struct(stat1,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do all manual steps for each mouse first:
if findUsvs
    for k = 1:totalMice
        mouseNum = mouseNums{1,k};
        waveFile = [wavPath, mouseNum, '.wav'];
        disp(['Reading WAV file: ', waveFile])
        [y, Fs] = audioread(waveFile);
        Yinfo= audioinfo(waveFile);
        totalSec =Yinfo.Duration;
        Sec90Samples = Fs*90;
        
        stopInd = length(y);
        startInd = stopInd - Sec90Samples;
        yTrim = y(startInd:stopInd); 
        calcDuration = (stopInd-startInd)/Fs;

        disp(['Taking FFT...'])
        nfft = 512; window = 512; noverlap = window*0.5;
        thresh = -95; % threshold in decibels;
        [~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs,'MinThreshold', thresh);  %do not use threshold if zscoring below
        %note that P is the power spectral density in W/Hz
        Tperiod =calcDuration/length(T);

        refPower = 10^-12; %eference power is 10^-12 watts (W), which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
        signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 
        signal(signal == -Inf) = 0; %set -Inf to zero for mean calculation below
        
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
        usvPowerPerSample = mean(signalCleaned); %mean power across frequencies (so total dB/Hz from lowFreq-highFreq, over temporal smoothing filter); do NOT normalize for now
        % for filtering/smooothing:
        wndSec = 0.1; %single USVs are on order of 50-100ms
        wndsz = round(wndSec/Tperiod); %convert seconds to samples
        gaussFilter = gausswin(wndsz); % wndsz value found by guess & check
        gaussFilter = gaussFilter ./ sum(sum(gaussFilter)); % normalize so that overall levels don't change
        validSamples = length(usvPowerPerSample)-wndsz+1;
         usvPowerPerSampleSmooth = conv(usvPowerPerSample, gaussFilter, 'valid'); %smoothing filter 

         %count USVs based on power peaks above threshold, which is somewhat arbitrary
         powerThresh = 2.3; 
        [pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
        numUsvs = length(locUsvs);
        usvPwrSum = sum(usvPowerPerSample);

%         figure;
%         subplot(2,1,1)
%         title('original');
%         imagesc(T, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
%         set(gca,'YDir','normal') %flip to make small freq on bottom
%         colormap(1-bone);
%         ylabel('frequency (kHz)');
%         xlabel('time (seconds)')
%         subplot(2,1,2)
%         title('power per sample');
%         plot(T(1:validSamples), usvPowerPerSampleSmooth)
%         hold on;
%         plot([0,calcDuration], [powerThresh, powerThresh], 'r--')
%         plot(locUsvs*Tperiod, ones(length(locUsvs),1).*max(pks), 'g*')
%         ylabel('USV power (mean dB/Hz in 40-100kHz band)');
%         xlabel('time (seconds)')
%         hold off;

        matName = [wavPath, mouseNum, '_usvData.mat']; %save data for each individual mouse
        save(matName, 'numUsvs', 'usvPwrSum');
    end
end

%% combine data across mice:
if combineMice
    for k = 1:totalMice
        mouseNum = mouseNums{1,k};
        matName = [wavPath, mouseNum, '_usvData.mat'];
        load(matName);
        usvStats.numUsvs(end+1) = numUsvs; 
    end
    
    save(matNameSaveSummary, 'usvStats');
end

%% Compute stats
if computeStats
    numGroupsToCompare = size(groupsToCompare, 2);

    for i = 1:numGroupsToCompare
        load([wavPath, groupsToCompare{1,i}, '.mat'])
         usvStatsVec(i,:) = usvStats.numUsvs;  %create matrix for comparisons
    end
    
     if normalize
        maxUsvs = max(usvStatsVec,[],1); %get individual max across days for each mouse
        for k=1:totalMice
            usvStatsVec(:,k) = usvStatsVec(:,k)./maxUsvs(k); % calculate to normalize by entire multi-day max
        end
     end
     
    [meanUsvs, stdUsvs, semUsvs] = grpstats(usvStatsVec',[],{'mean','std','sem'});

% plot data:
    cmap = colormap('lines');
    close(1); 
    fontSz = 20;
    
    %     plot mean # UVSs over days, with overlaid individual traces
    hOverDays = figure;
    hold on;
    for ip = 1:size(usvStatsVec,2)  %first plot individual mouse traces
        indTrace = [usvStatsVec(1,ip) usvStatsVec(2,ip) usvStatsVec(3,ip) usvStatsVec(4,ip) usvStatsVec(5,ip)  usvStatsVec(6,ip)];
        plot(indTrace, 'Color', cmap(ip,:), 'LineWidth', 2, 'LineStyle', '--')
    end
    errorbar([1 2 3 4 5 6], [meanUsvs(1) meanUsvs(2) meanUsvs(3) meanUsvs(4) meanUsvs(5) meanUsvs(6)], [semUsvs(1) semUsvs(2) semUsvs(3) semUsvs(4) semUsvs(5) semUsvs(6)], 'Color', 'k', 'LineWidth', 4);
    plot([1 2 3 4 5 6], [meanUsvs(1) meanUsvs(2) meanUsvs(3) meanUsvs(4) meanUsvs(5) meanUsvs(6)], 'o-', 'Color', 'k', 'MarkerSize', 5, 'LineWidth', 3);
    hold off;
    if normalize
        axis([0.5 6.5 0 1]);
    else
        axis([0.5 6.5 0 250]);
    end
    set(gca, 'FontSize', fontSz);
    set(gca,'XTick', [1 2 3 4 5 6], 'XTickLabel',{'p7';'p8';'p9';'p10';'p11';'p12'})
    ylabel('# USVs', 'FontSize', fontSz)
    %[p,tbl,stats] = kruskalwallis(x); %x is Sample data for the hypothesis test, specified as a vector or an m-by-n matrix. If x is an m-by-n matrix, each of the n columns represents an independent sample containing m mutually independent observations.
%     [pKruskalWallis,tblKruskalWallis,statsKruskalWallis] = kruskalwallis(usvStatsVec'); %try parameterized stats
%     cKruskalWallis = multcompare(statsKruskalWallis, 'estimate', 'kruskalwallis')
    %ksResult = kstest(usvStatsVec(1,:)) %test normality, returns 0 if normally distributed (repeat across mice, for each day)
    [pFriedman,tblFriedman,statsFriedman] = friedman(usvStatsVec',1); %nonparametric, based on ranks
    cFriedman = multcompare(statsFriedman, 'estimate', 'friedman', 'ctype', 'dunn-sidak')
end
