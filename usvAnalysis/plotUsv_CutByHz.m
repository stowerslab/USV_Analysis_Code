% top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior
% NOTE: requires Statistics toolbox
%     
%%
clear all; close all;
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-9 ChR2 female and controls\2018-9-10 ChR2 control and female_T1\';
%wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-9 ChR2 female and controls\2018-9-11 ChR2 control and female_T2\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-9 ChR2 female and controls\2018-9-13 CHR2 female_CONTROL_T4\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\All ChR2 wav files-FEMALE AND MALE\Male only\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-11 BNST Esr ChR2 batch2\2018-11-19  BNST ESR T6\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-11 BNST Esr ChR2 batch2\2018-11-16 BNST ESR CHR2_T4\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-11 BNST Esr ChR2 batch2\2018-11-14 BNST ESR CHR2_T2\';
% wavPath = 'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2018-11 BNST Esr ChR2 batch2\2018-11-20 BNST ESR T7\';
% wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior analysis\2019-7 GAQ vGat ChR2 hM3q fig and anlysis\calling animals_GAQ2 GAQ4 GAQ6\T1T2\';
% wavPath =  '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior data\CHR2 experiment\2019-6 EC Esr ChR2 Data-batch4\2019-6-10 EC T5 terminal time dose\';
% wavPath =  '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior analysis\2019-8 WT SMUF for USV paper\T1\';
% wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\CHR2 experiment\2019-10 PAG vGat ChR2\2019-10-27 PAG VGAT CHR2_T5\';
% wavPath =  '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior data\CHR2 experiment\2018-11TO2019-3 BNST Esr ChR2 batch2\2018-11-18 BNST ESR T5\';
%  wavPath =  '\\172.29.164.29\jingyi chen\Dropbox\Jingyi Data backup\Behavior data\CHR2 experiment\2019-10 VGAT CHR2 DREADD CONTROLS\2019-10-16 chr2_DREADD_control T2\';
%  wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\CHR2 experiment\2020-6 vGluT2 ChR2 LPOA\2020-5-22 vGluT2-ChR2-T2\';
% wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\DREADD experiment\2020-8 Esr DREADD control\2020-8-14 CNO day5\';
 wavPath =  'E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\Behavior data\2020-6 chr2 PAG-FP\2020-7-2 ChR2 T4\';
% vSheet = [wavPath '1-5-10HZ_Frame'];
% vSheet = [wavPath '25HZ_Frame'];
% vSheet = [wavPath '10HZ_Frame'];
% vSheet = [wavPath '5HZ_Frame'];
vSheet = [wavPath 'Timesheet_T4'];
% vSheet = [wavPath 'Time sheet_GAC_T1_B'];
stimLength = 5; %in seconds
fontSz = 25;
preSecToPlot = 10;
postSecToPlot = 20;%match GCaMP data, 10s before and after stimulation
% preSecToPlot = 5;
% postSecToPlot = 25;
% preSecToPlot = 0;% cut for DREADD data
% postSecToPlot = 240;% cut for DREADD data for 4 mins long
usvDataBndpwr = [];
    
plotSeparate = 0;

%read in stimulation timing data if necessary
vData = xlsread(vSheet);
videoRate = 15;
% videoRate = 14.997; %batch 2 has different video rate
% videoRate = 14.989 %For EC24_T6
% videoRate = 14.954 %For EC22_T4
% videoRate = 14.968 %For FEC19_T7
%videoRate = 14.998; %For EC18-24, FEC8-9 for T3

% mouseNums = {'EC8_T2','EC9_T2','EC10_T2','EC11_T2','EC12A_T2','EC15_T2',};
% mouseNums = {'FEC3_T1','FEC4_T1','FEC6_T1','FEC7_T1'};
% mouseNums = {'FEC3_T2','FEC4_T2','FEC6_T2','FEC7_T2'};
% % mouseNums = {'EC17_T2','EC18_T2','EC19_T2','EC20_T2','EC21_T2','EC22_T2','EC23_T2','EC24_T2','FEC8_T2','FEC9_T2'};
% % mouseNums = {'EC17_T4','EC18_T4','EC19_T4','EC20_T4','EC21_T4','EC22_T4','EC23_T4','EC24_T4','FEC8_T4','FEC9_T4'};
% % mouseNums = {'EC17_T5','EC18_T5','EC19_T5','EC20_T5','EC21_T5','EC22_T5','EC23_T5','EC24_T5','FEC8_T5','FEC9_T5'};
% % mouseNums = {'EC18_T6','EC19_T6','EC20_T6','EC21_T6','EC22_T6','EC23_T6','EC24_T6','FEC8_T6','FEC9_T6'};
% % mouseNums = {'EC17_T7','EC18_T7','EC19_T7','EC20_T7','EC21_T7','EC22_T7','EC23_T7','EC24_T7','FEC8_T7','FEC9_T7'};
% % mouseNums = {'EC24_T2'}; % do this one seperate because videorate is different for some animals
%mouseNums = {'EC18_T3', 'EC19_T3', 'EC20_T3', 'EC21_T3', 'EC22_T3', 'EC23_T3', 'EC24_T3', 'FEC8_T3', 'FEC9_T3'};
% mouseNums = {'FEC3_T3', 'FEC4_T3', 'FEC6_T3', 'FEC7_T3'};
% mouseNums = {'ESC1_T4', 'ESC2_T4', 'ESC3_T4', 'ESC4_T4' , 'ESC5_T4'};
% mouseNums = {'FEC10_T4','FEC11_T4','FEC12_T4'} ;
% mouseNums = {'FEC10_T5','FEC11_T5','FEC12_T5'} ;
% mouseNums = {'FEC10_T3UNI','FEC11_T3UNI','FEC12_T3UNI','GAD3_T3UNI'};
% mouseNums = {'FEC10_T3','FEC11_T3','FEC12_T3','GAD1_T3','GAD2_T3','GAD3_T3'};
% mouseNums = {'GAD1_S2_1-5-10HZ','GAD2_S2_1-5-10HZ','GAD3_S2_1-5-10HZ'};
% mouseNums = {'GAD1_S2_25HZ','GAD2_S2_25HZ','GAD3_S2_25HZ'};
% mouseNums = {'GAD1_S4_25HZ','GAD2_S4_25HZ','GAD3_S4_25HZ'};
% mouseNums = {'GAQ2_S4_10HZ','GAQ4_S4_10HZ','GAQ6_S4_10HZ'};
% mouseNums = {'GAQ4_S4_50HZ'};
% mouseNums = {'GAQ2_S4_50HZ20S','GAQ4_S4_50HZ20S','GAQ6_S4_50HZ20S'};
% mouseNums = {'GAQ2_T1_B','GAQ4_T1_B','GAQ6_T1_B'};
% mouseNums = {'GUC5_T3', 'GUC6_T3'};
% mouseNums = {'EC25_T3_B','EC26_T3_B','EC27_T3_B','EC28_T3_B'};
% mouseNums = {'EC25_T5','EC26_T5','EC27_T5','EC28_T5'};
% mouseNums = {'h1_1','h02_1','h2_1','h03_1','h3_1','h04_1','h4_1','M2_T1','M4_T1','M5_T1','M6_T1','M7_T1'}; % WT SMUF USV data
% mouseNums = {'GTCR3_T5'};
% mouseNums = {'WT1','WT2','WT3','WT4','WT5','WT6','WT7','WT8','WT9','WT10','WT11','WT12','WT13','WT14','WT15','WT16','WT17','WT18','WT19','WT20'};
% mouseNums = {'WT1_T2','WT2_T2','WT3_T2','WT4_T2','WT5_T2','WT6_T2','WT7_T2','WT8_T2','WT9_T2','WT10_T2','WT11_T2','WT12_T2','WT13_T2','WT14_T2','WT15_T2','WT16_T2','WT17_T2','WT18_T2','WT19_T2','WT20_T2'};
% mouseNums = {'GAC1_T2', 'GAC2_T2'};
% mouseNums = {'ED1_T1', 'ED2_T1','ED3_T1', 'ED4_T1','ED5_T1', 'ED6_T1','ED7_T1', 'ED8_T1'};
% mouseNums = {'ED1_C5', 'ED2_C5','ED3_C5', 'ED4_C5','ED5_C5', 'ED6_C5','ED7_C5', 'ED8_C5'};
% mouseNums = {'ED1_S4', 'ED2_S4','ED3_S4', 'ED4_S4','ED5_S4', 'ED6_S4','ED7_S4', 'ED8_S4'};
%  mouseNums = {'LUC1_T2', 'LUC2_T2','LUC3_T2', 'LUC4_T2'};
%   mouseNums = {'EDC1_C5', 'EDC2_C5','EDC3_C5', 'EDC4_C5','EDC5_C5'};
%    mouseNums = {'EDC1_S4', 'EDC2_S4','EDC3_S4', 'EDC4_S4','EDC5_S4'};
  mouseNums = {'GCG1_T4', 'GCG2_T4','GCG3_T4', 'GCG4_T4'};
% mouseNums = {'GAC3_T2','GAQ7_T2'};
vAllOffset = 0; %offset into video parameter spreadsheet45; for syncing to video only
totalMice = size(mouseNums, 2);

for j = 1:totalMice
    mouseNum = mouseNums{1,j};
    
    waveFile = [wavPath, mouseNum, '.wav'];
    disp(['Reading WAV file: ', waveFile])
    [y, Fs] = audioread(waveFile);
    i = audioinfo(waveFile); % i.TotalSamples, i.Duration, i.SampleRate

    syncVidSample = vData(vAllOffset+j,1);
    stimStartVidSamples = vData(vAllOffset+j,2:end) - syncVidSample; %subtracted sync frame from stim frames
    stimStartVidSec = stimStartVidSamples./videoRate;
    
    preAudioSamplesToPlot = preSecToPlot*i.SampleRate;
    postAudioSamplesToPlot = postSecToPlot*i.SampleRate;
    stimStartAudioSamples = stimStartVidSec*i.SampleRate; %assumes file starting at red beep
        %% choose threshold for USVs in frequency domain
    display('Getting USV threshold...')
    % Fill in header information (fake it where necessary):******************************
    nfft = 4096; %sets frequency resolution; use higher res for finding sync pulse, but this is reset below for vocalizations
    options.numoverlap = nfft/2;  %note that this will set temporal resolution ~(total samples)/(nfft-overlap), so smaller nfft & larger overlap means more temporal resolution
    header.nscans = i.TotalSamples;
    header.duration = i.Duration; %length of recording in seconds
    header.scanrate = i.SampleRate;
    header.tacq = header.nscans/header.scanrate;  %this is also the number of seconds in recording
    header.numch = i.NumChannels;
    header.channels = 0:i.NumChannels-1;
    header.scalemult = 1;
    header.scaleoff = 0;
    header.voltageMin = -Inf;
    header.voltageMax = Inf;
    header.date = '';
    header.time = '';
    header.usrhdr = '';
    header.flag = 'WU1';    
    header.nfreq = nfft/2+1;
    header.memmax = 5*1024*1024;
    header.nblocks = floor((header.nscans-options.numoverlap)/(nfft - options.numoverlap));
    header.nptsPerBlock = header.nscans/header.nblocks;
    header.columnTotal = (header.nptsPerBlock/header.nfreq)*header.nblocks; % No "2" because half-overlap
    % header.freqMin = freqRange(1);
    % header.freqMax = freqRange(2);
    % **********************************************************************************

    % usvThresh = chooseUsvThresh(y,Fs,nfft) %maunally find threshold (in freq. domain) for USVs first time with new recording setup
    usvThresh = 1;
    header.threshold = usvThresh;

    %% sync to video using Holy code
    display('Syncing to video...')
    thresh = 20;  % threshold in freq. domain for use with sync pulse; found by trial and error changed by Jingyi to 18 on 6/7/19 used to be 7 and can't detect accurately for high noise audios
    %use 3 for the new beep since it is much quieter
    options.band = [4 5]; %in KHz, for finding sync pulse (s/b ~4.6kHz);
    secLookForSync = 10; %cut first 10 seconds only to find sync pulse
    cutTime = header.scanrate*secLookForSync;  
    ySync = y(1:cutTime);
    [sng,f,t] = sparsesng(ySync,Fs,nfft,thresh,options); %#ok<ASGLU> %use Holy sparse representation

    % from spsngplot: first take complex FFT values are convert to real lognormal values for plotting
    [~,j,~] = find(sng);  %find nonzero elements of the sparse matrix, put in 's'
    beepStart = t(min(j)); %j is the index into t at which the first above-thresh sample occurs, which is the sync pulse
%     beepStart = 5.5; % change it here manually if can't use code to find the start
    if isempty(beepStart)
        beepStart =  4; %if no beep, just assume 4 sec offset
        display('WARNING!!!!!!!!!!!: Red beep not found; assumning 2 sec...')
    end

%     % ALTERNATIVE manual selection of first sync sample
%     bandsig = find(f >= options.band(1) & f <= options.band(2)); %subset of freqs of interest
%     fPlot = linspace(options.band(1),options.band(2),length(bandsig));
%     hFig = figure;
%     spsngplot(hFig,sng,fPlot,bandsig,t,' Zoom in, then click on start of sync tone:');
%     
%     % now manually select time point corresponding to start of red beep:
%     if v.startRedBeep ~= 0  %if there is a red beep, find it on spectrogram
%         hAxes = gca;  % get handle to axes for 'waitfor' funtion below
%         zoom on
%         waitfor(hAxes,'CameraPosition');  % wait until zoom changes - note that you only have one chance to change ;-)
%         zoom off
%         hBeepStart = impoint(hAxes);  %first draw ellipse center point manually
%         pos = getPosition(hBeepStart);
%         beepStartMin = pos(1);  %start time of beep in seconds
%     else
%         beepStartMin =  4; %if no beep, just assume 4 sec offset
%     end
%     close(hFig); %close sync finder figure
%     pause on; pause(0.1); pause off;

    %% calculate range of data within assay, using video parameters
    [~, cutStartBlock] = min(abs(t - beepStart)); %find closest FFT index to time point chosen
    cutStartSample = floor(cutStartBlock*header.nptsPerBlock);
    % cutStopSample = round(length(y)*0.5);
    % keyboard
    yTrim = y(cutStartSample:end);
    clear y;

    %% take FFT & plot
    display('Taking FFTs & plotting spectra...')
    nfft = 2048;
    window = 1024; %in samples at 192KHz (1000 ~= 5ms)
    noverlap = window*0.5;
    thresh = -90; %in decibels; this is included in spectrogram 'MinThreshold' in newer versions Jingyi changed to -90 instead of -80 on 7/9/2019
    lowFreq = (45/96)*nfft/2; %'xx/96' in kHz
    highFreq = (96/96)*nfft/2; % CHANGE FROM 80 to 96/96 by Jingyi 2019-2-14
    % figure;
%      for i = [1]
%     for i = [1 2 3 4 5 6 7 8 9 10]
    for i = [1 2 3 4 5 6 7 8 9 10 11 12 13] %choose the columns of interest; 1:length(stimStartAudioSamples)
% for i = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]
%     for i = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20] %choose the columns of interest; 1:length(stimStartAudioSamples), for T8 and T9 of unilateral
        yCurrent = yTrim((stimStartAudioSamples(i)-preAudioSamplesToPlot):(stimStartAudioSamples(i)+postAudioSamplesToPlot));
        [~,F,T,P] = spectrogram(yCurrent,window,noverlap,nfft,Fs);
        
        filname= [mouseNum, '_column', int2str(i), '.wav'];
        audiowrite(filname,yCurrent,Fs)
%         
        signal = 10*log10(abs(P));
        signal(signal<thresh) = thresh;
        signalBndpwr = mean(signal(lowFreq:highFreq,:),1);% ./ max(signal(:,:),[],1);
        usvDataBndpwr(:,end+1) = signalBndpwr - signalBndpwr(1); %subtract starting value
        timescale = T - preSecToPlot;
        x1 = [timescale(find(timescale>0, 1, 'first')) timescale(find(timescale>0, 1, 'first'))];
        x2 = [timescale(find(timescale>stimLength, 1, 'first')) timescale(find(timescale>stimLength, 1, 'first'))];

        if plotSeparate
            figure; %subplot(length(stimStartAudioSamples),1,i);
            imagesc(timescale, F(lowFreq:highFreq)./1000, signal(lowFreq:highFreq,:)); 
            set(gca,'YDir','normal') %flip to make small freq on bottom
            colormap(1-gray);
%             caxis([-thresh -50]);
            ylabel('frequency (kHz)');
            xlabel('time (seconds)')
            set(gca,'XTick',[x1(1) x2(1)])
            set(gca,'XTickLabel',{'0' num2str(stimLength)})
            yLims = get(gca, 'YLim');
            y1 = [yLims(1) yLims(2)];
            hold on;
            plot(x1, y1, 'b--', 'LineWidth', 2); %line for stimulus start
            plot(x2, y1, 'b--', 'LineWidth', 2); %line for stimulus end
            set(gca, 'FontSize', 20);
            hold off;
        end
    end

end

% hBandPowerHeatmap = figure;
% % [~,maxRowIndeces] = sort(mean(usvDataBndpwr)); %sort by max mean bandpower
% % usvDataBndpwrSorted = usvDataBndpwr(:,maxRowIndeces);
% x1samp = [find(timescale>0, 1, 'first') find(timescale>0, 1, 'first')];
% x2samp = [find(timescale>stimLength, 1, 'first') find(timescale>stimLength, 1, 'first')];
% set(gca,'XTick',[x1samp(1) x2samp(1)])
% set(gca,'XTickLabel',{'0' num2str(stimLength)})
% heatmapOld(usvDataBndpwr')
% yLims = get(gca, 'YLim');
% y1 = [yLims(1) yLims(2)];
% hold on;
% plot(x1samp, y1, 'b--', 'LineWidth', 2); %line for stimulus start
% plot(x2samp, y1, 'b--', 'LineWidth', 2); %line for stimulus end
% ylabel('mean power spectral density 45-85kHz (dB/Hz)', 'FontSize', fontSz);
% xlabel('time (sec)', 'FontSize', fontSz);
% caxis([0 1]);
% colormap(1-gray)
% colorbar
% set(gca, 'FontSize', fontSz);
% %hold off;
