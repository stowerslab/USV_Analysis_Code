% top-level script to process USV audio data for scent-mark-to-female-urine (SMUF) behavior
% NOTE: requires Statistics toolbox

%%
clear all; close all;
wavPath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior data\2018-9 GCaMP BNST\2018-9-24 GCaMP T5\';

vSheet = [wavPath 'Trim sheet_T5'];
fontSz = 20;
usvDataBndpwr = [];
    
%read in stimulation timing data if necessary
vData = xlsread(vSheet);
% videoRate = 14.499;  %Have to change based on AE's reads on video rate


%  mouseNums = {'BG9_T1' , 'BG10_T1', 'BG11_T1', 'BG12_T1', 'BG13_T1', 'BG14_T1'};
mouseNums = {'BG6_T5L' ,'BG7_T5L' ,'BG8_T5L'};
% mouseNums = {'BG6_T6R' ,'BG7_T6R' ,'BG8_T6R'};
vAllOffset = 0; %offset into video parameter spreadsheet; for syncing to video only
totalMice = size(mouseNums, 2);

for j = 1:totalMice
    mouseNum = mouseNums{1,j};
    
    waveFile = [wavPath, mouseNum, '.wav'];
    display(['Reading WAV file: ', waveFile])
    [y, Fs] = audioread(waveFile);
    i = audioinfo(waveFile); % i.TotalSamples, i.Duration, i.SampleRate
    videoRate = vData(vAllOffset+j,6); % read video rate, might need to change depend on excel
    syncVidSample = vData(vAllOffset+j,1);
    stimStartVidSamples = vData(vAllOffset+j,2) - syncVidSample; %subtracted sync frame from start frames
    stimStartVidSec = stimStartVidSamples./videoRate;
    
    stimEndVidSamples = vData(vAllOffset+j,3) - syncVidSample; %subtracted sync frame from end frames
    stimEndVidSec = stimEndVidSamples./videoRate;
    
    
    stimStartAudioSamples = stimStartVidSec*i.SampleRate; %assumes file starting at red beep
    stimEndAudioSamples = stimEndVidSec*i.SampleRate;
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
    thresh = 7;  % threshold in freq. domain for use with sync pulse; found by trial and error
    options.band = [4 5]; %in KHz, for finding sync pulse (s/b ~4.6kHz);
    secLookForSync = 15; %cut first 15 seconds only to find sync pulse
    cutTime = header.scanrate*secLookForSync;  
    ySync = y(1:cutTime);
    [sng,f,t] = sparsesng(ySync,Fs,nfft,thresh,options); %#ok<ASGLU> %use Holy sparse representation

    % from spsngplot: first take complex FFT values are convert to real lognormal values for plotting
    [~,j,~] = find(sng);  %find nonzero elements of the sparse matrix, put in 's'
    beepStart = t(min(j)); %j is the index into t at which the first above-thresh sample occurs, which is the sync pulse
    if isempty(beepStart)
        beepStart =  4; %if no beep, just assume 4 sec offset
        display('WARNING!!!!!!!!!!!: Red beep not found; assumning 4 sec...')
    end

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
    thresh = -80; %in decibels; this is included in spectrogram 'MinThreshold' in newer versions
    lowFreq = (45/96)*nfft/2; %'xx/96' in kHz
    highFreq = (85/96)*nfft/2;
    % figure;

%     for i = [1 2] %choose the columns of interest; 1:length(stimStartAudioSamples)
        yCurrent = yTrim(stimStartAudioSamples:stimEndAudioSamples);
        [~,F,T,P] = spectrogram(yCurrent,window,noverlap,nfft,Fs);
        
        filname= [mouseNum, '_FP', '.wav'];
        audiowrite(filname,yCurrent,Fs)
%     end
end
