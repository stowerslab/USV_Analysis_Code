%% Header
clear all; close all;

% need to be in directory with data E:\Jingyi\Dropbox (Scripps Research)\Jingyi Data backup\GCaMP Data\BNST master gcamp
% this makes the averaged gcamp signal at sentence timepoints with a
% scrambled baseline gcamp sig for comparison
section_type = 'with Female'; % baseline, spontaneous, with Female, or with FU; if no syllables might error

SAVE_SIGNAL_PATH = 'E:\Pete\2020-7 PAG vglut2 gcamp + audio data\'; % where to save the signal variables after processing is done
sync_needed = 1;

wavPath = [pwd '\' section_type ' wav\'];

gcampPath = [pwd '\' section_type '\'];
% no_usv_gcampPath = [pwd '\baseline\'];

% for syncing wav and gcamp need a timings excel sheet
if sync_needed == 1
    timings = readtable('timings.xlsx');
    firstCol = cell2mat(timings(:,1)); % trial names col from table
end

WAVFILES = dir(fullfile(wavPath,'*.wav')); %gets all wav files in struct

GCAMPFILES = dir(fullfile(gcampPath,'*.mat')); 

% noUSV_GCAMPFILES = dir(fullfile(no_usv_gcampPath, '*.mat'));

% some variables for later:
cutoff_df = 1.5; % if max of gcamp signal is less than cutoff, exclude it from being processed

skipped_IDs = []; % matrix to hold id's of skipped trials, write to txt file at end

vidFPS = 15; % for calculating sync diff between audio and gcamp

all_bout_lengths = []; % store bout lengths or sentence length. average last locUSV - first locUSV in a sentence
init_bout_lengths = []; 
first3_bout_lengths = [];
saved_IDs = [];

gcampFrameRate = 20;
%extended time (s) from onset to draw figure
timeExtension = 10;
%a local time (s) window from onset to search for peaks
localWindow = 0.25;

minSyllableDistance = 2; %define a sentence by 2s breaks

%store processed data
avg_sig_norm_cleaned = []; % zeros(length(GCAMPFILES),length(-timeExtension:1/gcampFrameRate:timeExtension));
avg_scrambledSig_norm = []; % zeros(length(GCAMPFILES),length(-timeExtension:1/gcampFrameRate:timeExtension));
no_usv_sig_norm = [];
no_usv_Same_Section = []; % pick random timepoints from same section that aren't usv timepoints

k = 0; % just a counter/indexer in the for loop, used to index into gcampCell cell-array


for wavFile = WAVFILES' % loop through wav files, match to its gcamp file, combine all in the end
    
    % get id from wav file
    
    mouseID = strrep(wavFile.name, '-', '_');
    
%     expression = '[A-Z]{2}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}'; % regex expression to match to the mouseid in the wav filename
    expression = '[A-Z]{2}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}_[A-Z]{1}'; % PAG id regex
%     expression = '[A-Z]{3}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}_[A-Z]{1}'; % POA id regex

    mouseID = regexp(mouseID,expression,'match');  % wav mouse ID

    for gcampFile = GCAMPFILES'
        
        % get id's from gcamp files
        gMouseID = strrep(gcampFile.name, '-', '_');
    
%         expression = '[A-Z]{2}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}'; % regex expression to match to the mouseid in the wav filename
        expression = '[A-Z]{2}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}_[A-Z]{1}'; % PAG id regex
%         expression = '[A-Z]{3}[0-9]{1,2}_[A-Z]{1}[0-9]{1,2}_[A-Z]{1}'; % POA id regex
        gMouseID = regexp(gMouseID,expression,'match'); % gcamp mouse ID
        
        if strcmp(gMouseID,mouseID)
            load(fullfile(gcampPath,gcampFile.name));
            
            % TRIM USING TIMINGS SHEET
            if sync_needed == 1
                row = timings(strcmp(firstCol, gMouseID), :);
                if isempty(row)

                    continue;
                end
                
                % 1. frame difference between the two signals from the 15fps video
                frameDiff = row.GCAMPFrame - row.USVWAVFrame ; 
                secDiff = frameDiff/vidFPS; % USE THIS TO TRIM AUDIO FROM BEGINNING

                % 2. if the timings in excel are seconds then use this
%                 secDiff = row.GCAMPTime - row.USVWAVTime; 
            end
            
            disp(['Reading GCaMP file: ', gcampFile.name]);
            sig_norm = sig_norm * 100;
            
            skipped = 0;
            if max(sig_norm) < cutoff_df % checks if the entire gcamp signal has a value higher than cutoff
                disp(['SKIPPED max df: ' num2str(max(sig_norm))]);
                skipped_IDs = [skipped_IDs; gMouseID];
                continue;
            else
                sig_norm = zscore(sig_norm);
                
                k = k + 1;
                %% USV Section

                AUCpowerARRAY = [];
                usvPwrPerSampleAll = [];    

                %Read USV file that has been trimmed to match GCamP recording time

                disp(['Reading WAV file: ', wavFile.name]);

                [y, Fs] = audioread(fullfile(wavPath,wavFile.name));
                i = audioinfo(fullfile(wavPath,wavFile.name));

                totalSec = i.Duration;
                USVSamplerate=i.SampleRate;

%                 startInd = 1;
                % ASSUMING AUDIO STARTS BEFORE GCAMP
                if sync_needed == 1
                    startInd = secDiff*USVSamplerate; % where the sync start should be
                else
                    startInd = 1;
                end
                stopInd = startInd + length(sig_norm)/gcampFrameRate * USVSamplerate; % to match time, should be less than 0.5 s different

                if stopInd > size(y,1) % if the audio file is shorter than gcamp, trim gcamp to audio length
                    stopInd = size(y,1);
                    sig_norm = sig_norm(1:round((size(y,1)-startInd)/USVSamplerate * gcampFrameRate));
                end
                
                yTrim = y(startInd:round(stopInd)); % trim audio to gcamp length and synced start
                clearvars y;
                
                %%

                GCaMPCell{k}=sig_norm;
                numGcampSamples = length(GCaMPCell{k});
                totalGCampSec= length(sig_norm)/gcampFrameRate;
                timescaleGcamp= 0:1/gcampFrameRate:totalGCampSec-(1/gcampFrameRate);
            %     GCaMPtimesCell{k} = 0:gcampFrameRate:(numGcampSamples*gcampFrameRate)-gcampFrameRate; % 



                disp('Taking FFT...');
                nfft = 512;
                window = 512;
                noverlap = window*0.5;
                % thresh = -90; % threshold in decibels;
                [~,F,T,P] = spectrogram(yTrim,window,noverlap,nfft,Fs); %,'MinThreshold', thresh);  %do not use threshold if zscoring below
                %note that P is the power spectral density in W/Hz
                Tperiod =i.Duration/length(T);
                

                refPower = 10^-12; %reference power is 10^-12 watts (W)/m^2, which is the lowest sound persons of excellent hearing can discern (wiki, http://www.sengpielaudio.com/calculator-soundpower.htm)
                signal = 10*log10(abs(P./refPower)); %convert to dB for acoustic convention (now signal is dB/Hz) 
                clearvars T P;
                disp('Filtering noise...')
                % idea here take z-score & compare to other freqs to remove broadband noise: 
                zsignal = zscore(signal);
                lowFreq = find(F>40000,1,'first');  % index for lowpass cutoff freq
                highFreq = find(F>90000,1,'first'); % index for high cutoff used below --- jingyi says 90khz, og: 80khz - pete
                zsignal(1:lowFreq,:) = 0; % lowpass - set everything below cutoff to zero
                zsignal(highFreq:end,:) = 0; % highpass - add by Jingyi on 11/29/2018
                zthresh = 1.5; %1.3 for pup USV
                zsignal(zsignal<zthresh) = 0; %threshold zscore
                signalCleaned = signal; %create a copy to clean below
                clearvars signal;
                signalCleaned(zsignal==0) = 0; %JAK find where zscore reduced noise and artificially set that back into original file (so unit still dB/Hz); could use morphological expansion here to be more conservative!!!
                clearvars zsignal;
                
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
                samplerateRatio= length(usvPowerPerSampleSmooth)/length(sig_norm);
                timescaleUSV=0:1/(gcampFrameRate*samplerateRatio):totalGCampSec-1/(gcampFrameRate*samplerateRatio); % match USV time with GCamP time

                AUCpower(k) = trapz(usvPowerPerSample)/totalSec; 
                powerThresh = 1; % 
%                 figure;
%                 findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150)
                [pks, locUsvs] = findpeaks(usvPowerPerSampleSmooth,'MinPeakHeight', powerThresh, 'MinPeakDistance', 0.150); %peaks must be separated by ~150ms
                numUsvs(k) = length(locUsvs); %number of USVs

                %find syllables with minimum distance
                syllable = zeros(length(timescaleUSV),1);
                clearvars timescaleUSV;
                syllable(locUsvs,1)=1;
                
                % locSyllable is the sentence onset syllable
                
                [~,locSyllable] = findpeaks(syllable,'MinPeakDistance',minSyllableDistance/Tperiod);
                clearvars syllable;
                
                %% BOUTLENGTH
                % find start and end of sentences!
                if length(locSyllable) > 1
                    sentenceLocs=ismember(locUsvs, locSyllable); % logical array
                    sentenceLocs=double(sentenceLocs); % convert to numeric
                    for idx = find(sentenceLocs==1)
                        if idx > 1
                            sentenceLocs(idx-1) = 2; % the start locations have 1's, set the end locations to 2 
                        end
                    end
                    % we know locSyllable has the start locations, need to
                    % extract only the ends.
                    endLocs = transpose([locUsvs(find(sentenceLocs==2)), locUsvs(end)]); % locations of end of sentence

                    sentenceLength = endLocs - locSyllable;
                    all_bout_lengths = [all_bout_lengths; sentenceLength*Tperiod]; % append in and conert to seconds
                    init_bout_lengths = [init_bout_lengths, (endLocs(1) - locSyllable(1))*Tperiod];
                    if length(locSyllable) > 2
                        first3_bout_lengths = [first3_bout_lengths; (endLocs(1:3) - locSyllable(1:3))*Tperiod];
                    end
                end
                
                % time extension causing negative indexing, fixed by having minimum index = 1, and then lengthening the result array by zeroes to proper length - pete
                % ORIGINAL, USE ALL SYLLABLE ONSETS
                %                 all_onset_sig_norm = NaN(length(locSyllable),length(-timeExtension:1/gcampFrameRate:timeExtension));
                % MAY 19: USE INITIAL SYLLABLE ONLY
%                 locSyllable has onset timings for each syllable
                all_onset_sig_norm = NaN(1,length(-timeExtension:1/gcampFrameRate:timeExtension));
                
                %% Modify the for loop to either include gcamp at all syllables or only 1 syllable
                % TAKE ONLY THE FIRST SYLLABLE OF EACH PHRASE 5/19/2020 by
                % changing to "for n = 1:1" instead of "for
                % n=1:length(locSyllable)"

                for n = 1:1%length(locUsvs)  % <- the for-loop to change
                    
                    %preid is the closest time point in gcamp for each onset
                    if length(locUsvs) < n % for 5/19/2020
                        
                        break;
                    end
                    % timescaleGCAMP is an array, subtract all positions by
                    % the location of the syllable and get the abs min to
                    % get matching INDEX of syllable in timescalegcamp
                    
                    [~,preid] = min(abs(timescaleGcamp-locUsvs(n)*Tperiod));  % preid is a gcamp index
                    
                    %find the most prominent peak within a local time window to correct misalignment
                    
                    % cut gcamp
                    if preid+localWindow*gcampFrameRate < size(sig_norm, 2) 
    %                         disp(preid-localWindow*gcampFrameRate);
    %                         disp(preid+localWindow*gcampFrameRate);
                        [~,locOnset,~,prm] = findpeaks(sig_norm(1,max(1,preid-localWindow*gcampFrameRate):preid+localWindow*gcampFrameRate));
                        if length(locOnset) == 0
                            locOnset = [1];
                            prm=[1];
                        end
                    else
                        [~,locOnset,~,prm] = findpeaks(sig_norm(1,max(1,preid-localWindow*gcampFrameRate):length(sig_norm)));
                        if length(locOnset) == 0
                            locOnset = [1];
                            prm=[1];
                        end
                    end
                    
%                     if length(locOnset)>=1 && (preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 +timeExtension*gcampFrameRate) < size(sig_norm, 2) 
                     if length(locOnset)>=1 && (preid + locOnset(find(prm==max(prm))) -1) < size(sig_norm, 2) 
%                         correct misalignment
%                         postid = preid - localWindow*gcampFrameRate + locOnset(find(prm==max(prm))) -1 ;
                        
                        B = sig_norm(1, max(1,preid-timeExtension*gcampFrameRate):min(length(sig_norm),preid+timeExtension*gcampFrameRate));
                        
%                         if B(1) == 0
%                             disp(preid)
%                         end
%                         B = sig_norm(1, max(l,preid-timeExtension*gcampFrameRate):preid+timeExtension*gcampFrameRate);
                        % all_onset_sig_norm(n,:) = [zeros(1, size(all_onset_sig_norm, 2) - size(B, 2)), B]; % ISSUE: this moves signal to the right
                        
                        % pad differently if not using initial syll, since
                        % initial will always have nans appear before
                        if length(B) < length(all_onset_sig_norm) % IF TOO SHORT, pad with nan
                            if preid+timeExtension*gcampFrameRate > length(sig_norm) % padding END
                                
                                all_onset_sig_norm(n,:) = [B, NaN(1,preid+timeExtension*gcampFrameRate-length(sig_norm))];
                            else % padding BEGINNING
                                
                                beginIdx = -1*(preid-timeExtension*gcampFrameRate-1); % where the signal should start
    %                             disp('begin')
    %                             disp(length(B));
    %                             disp(beginIdx);
    %                             ix = floor(lendiff/2);
                                all_onset_sig_norm(n,:) = [NaN(1,round(beginIdx)), B];%length(all_onset_sig_norm) - length(B);
    %                           all_onset_sig_norm(n,floor(lendiff/2)+1:floor(lendiff/2)+length(B)) = B;
                            end
                            
                            
                        else % no need to pad with nans
                            
                            all_onset_sig_norm(n,:) = B;
                        end
%                         disp(all_onset_sig_norm(n,:));
%                         disp(preid);
                        
                        
                    end
                end
                
%                 disp(length(locSyllable));

                
            end

            saved_IDs = [saved_IDs; gMouseID];
%             avg_sig_norm_cleaned = [avg_sig_norm_cleaned;nanmean(all_onset_sig_norm(all_onset_sig_norm(:,1)~=0,:))];
            % MAY/19/2020 CHANGE BACK TO COMMENTED ABOVE TO GET
            % ORIGINAL. this for initial syllables only
            avg_sig_norm_cleaned = [avg_sig_norm_cleaned;(all_onset_sig_norm(all_onset_sig_norm(:,1)~=0,:))];
            
            
            %% PICK RANDOM NO-USV TIMEPOINTS FROM SAME SECTION
            num_no_usv = length((all_onset_sig_norm(all_onset_sig_norm(:,1)~=0,:))); % how many timepoints to pick from no-usv
            SS_random_NOUSV_id = datasample(round(1+timeExtension*gcampFrameRate):length(sig_norm(1,1:length(sig_norm)-timeExtension*gcampFrameRate)), num_no_usv);
            SS_random_NOUSV_id = setdiff(SS_random_NOUSV_id, locUsvs);
            if length(SS_random_NOUSV_id) < num_no_usv
                while 1
                    append_NOUSV_ids = datasample(round(1+timeExtension*gcampFrameRate):length(sig_norm(1,1:length(sig_norm)-timeExtension*gcampFrameRate)), num_no_usv-length(SS_random_NOUSV_id)  );
                    SS_random_NOUSV_id = [SS_random_NOUSV_id, append_NOUSV_ids];
                    SS_random_NOUSV_id = setdiff(SS_random_NOUSV_id, locUsvs);
                    if length(SS_random_NOUSV_id) == num_no_usv
                        break
                    end
                end
            end
            
            SS_nousv_scrambledSig_norm = nan(num_no_usv,length(-timeExtension:1/gcampFrameRate:timeExtension));
            SS_nousv_scrambledSig_norm = [nan(1, length(-timeExtension:1/gcampFrameRate:timeExtension)); SS_nousv_scrambledSig_norm];

            if num_no_usv < 1 % skip if no syllables
                continue;
            end

            for r = 1:length(SS_random_NOUSV_id)
                C = sig_norm(1,max(1,SS_random_NOUSV_id(r)-timeExtension*gcampFrameRate):SS_random_NOUSV_id(r)+timeExtension*gcampFrameRate);

                if length(C) < size(SS_nousv_scrambledSig_norm, 2)
                    continue;
                else
                    SS_nousv_scrambledSig_norm(r,:) = C;
                end
            end

            no_usv_Same_Section = [no_usv_Same_Section; SS_nousv_scrambledSig_norm];
            clearvars C;
            
            
            
            %% scramble timepoints
            numScrambled = 200; %size(all_onset_sig_norm(all_onset_sig_norm(:,1)~=0,:),1); % 
            scrambledID =  datasample(1:length(sig_norm(1,1:numGcampSamples-timeExtension*gcampFrameRate)),numScrambled);
            scrambledSig_norm = zeros(numScrambled,length(-timeExtension:1/gcampFrameRate:timeExtension));
            
            for j = 1:numScrambled
                %disp(max(1,scrambledID(j)-timeExtension*gcampFrameRate));
                %disp(scrambledID(j)+timeExtension*gcampFrameRate);

                % C = scrambled signal segment
                C = sig_norm(1,max(1,scrambledID(j)-timeExtension*gcampFrameRate):scrambledID(j)+timeExtension*gcampFrameRate);
                % pad the random selection with nan if too short

                if length(C) < size(scrambledSig_norm, 2)
%                             beginIdx = -1*(scrambledID(j)-timeExtension*gcampFrameRate-1);
%                             scrambledSig_norm(j,:) = [NaN(1,beginIdx), B];
%                             lendiff = length(scrambledSig_norm) - length(C);
%                             ix = floor(lendiff/2);
% %                             scrambledSig_norm(j,floor(lendiff/2)+1:floor(lendiff/2)+length(B)) = B;
%                             scrambledSig_norm(j,:) = [C,nan(1, size(scrambledSig_norm, 2) - size(C, 2))];
                    continue;
                else
                    scrambledSig_norm(j,:) = C;
                end

            end
            avg_scrambledSig_norm = [avg_scrambledSig_norm;nanmean(scrambledSig_norm)];
            clearvars C;
                % MAY/19/2020
%                         avg_scrambledSig_norm = [avg_scrambledSig_norm;(scrambledSig_norm)];
            
        else
%             skipped_IDs = [skipped_IDs; gMouseID];
            continue;
        end

        
    % load corresponding no-USV gcamp file
%         for noUSVgcampFile = noUSV_GCAMPFILES'
%             
%             gcamp_ID = split(gcampFile.name, 'processed');
%             gcamp_ID = gcamp_ID(1);
%             no_usv_gcamp_ID = split(noUSVgcampFile.name, 'processed');
%             no_usv_gcamp_ID = no_usv_gcamp_ID(1);
%             
%             
%             
%             if strcmp(gcamp_ID,no_usv_gcamp_ID)
%                 disp(noUSVgcampFile.name);
%                 load(fullfile(no_usv_gcampPath,noUSVgcampFile.name));
%                 sig_norm = sig_norm * 100;
%                 sig_norm = zscore(sig_norm);
%                 
%                 % pick random time points from corresponding control gcamp
%                 % signal and get gcamp signal at those times.
%                 num_no_usv = length((all_onset_sig_norm(all_onset_sig_norm(:,1)~=0,:))); % how many timepoints to pick from no-usv
%                 random_NOUSV_id = datasample(round(1+timeExtension*gcampFrameRate):length(sig_norm(1,1:length(sig_norm)-timeExtension*gcampFrameRate)), num_no_usv);
%                 
%                 nousv_scrambledSig_norm = nan(num_no_usv,length(-timeExtension:1/gcampFrameRate:timeExtension));
%                 nousv_scrambledSig_norm = [nan(1, length(-timeExtension:1/gcampFrameRate:timeExtension)); nousv_scrambledSig_norm];
%                 
%                 if num_no_usv < 1 % skip if no syllables
%                     continue;
%                 end
%                 
%                 for r = 1:length(random_NOUSV_id)
%                     C = sig_norm(1,max(1,random_NOUSV_id(r)-timeExtension*gcampFrameRate):random_NOUSV_id(r)+timeExtension*gcampFrameRate);
%                     
%                     if length(C) < size(nousv_scrambledSig_norm, 2)
%                         continue;
%                     else
%                         nousv_scrambledSig_norm(r,:) = C;
%                     end
%                 end
%                 
%                 no_usv_sig_norm = [no_usv_sig_norm; nousv_scrambledSig_norm];
%                 clearvars C;
%             end
%         end

    end

    
end
% 
% 
disp(['Saving signal variables: ' section_type]);
save([SAVE_SIGNAL_PATH section_type '.mat'],'timeExtension','gcampFrameRate', 'section_type', 'localWindow', 'avg_scrambledSig_norm', 'avg_sig_norm_cleaned');

disp('Plotting...');
figure;
% p0 = plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned',':k');
hold on;
title([section_type ' Z-Score GCaMP signals (initial syll)']);
% ylabel ('dF/F (%)');
ylabel ('Avg. Z-Scored dF/F');
xlabel ('time (s)');
% p1 = plot(-timeExtension:1/gcampFrameRate:timeExtension,mean(avg_sig_norm_cleaned, 1),'g');
p1 = boundedline(-timeExtension:1/gcampFrameRate:timeExtension,nanmean(avg_sig_norm_cleaned, 1), nanstd(avg_sig_norm_cleaned)/sqrt(size(avg_sig_norm_cleaned, 1)),'g'); % standard error of the mean is usually estimated as the sample standard deviation divided by the square root of the sample size
p2 = plot(-timeExtension:1/gcampFrameRate:timeExtension,nanmean(avg_scrambledSig_norm, 1),'k');
legend([p1, p2],'Average GCaMP','Average Scrambled GCaMP');

savefig([section_type ' ( ' num2str(cutoff_df) '% cutoff) ' num2str(timeExtension) 's no alignment std avg Z-score initial syll.fig']);
% savefig([section_type ' ( ' num2str(cutoff_df) '% cutoff) ' num2str(timeExtension) 's no alignment std avg initial syll.fig']);




% % ZSCORE SCALE
% figure;
% % p0 = plot(-timeExtension:1/gcampFrameRate:timeExtension,avg_sig_norm_cleaned',':k');
% hold on;
% title([section_type ' GCaMP signals (all syll)']);
% ylabel ('Z-score');
% xlabel ('time (s)');
% % p1 = plot(-timeExtension:1/gcampFrameRate:timeExtension,mean(avg_sig_norm_cleaned, 1),'g');
% p1 = boundedline(-timeExtension:1/gcampFrameRate:timeExtension,zscore(nanmean(avg_sig_norm_cleaned, 1)), nanstd(avg_sig_norm_cleaned)/sqrt(size(avg_sig_norm_cleaned, 1)),'g'); % standard error of the mean is usually estimated as the sample standard deviation divided by the square root of the sample size
% % p2 = plot(-timeExtension:1/gcampFrameRate:timeExtension,zscore(nanmean(avg_scrambledSig_norm, 1)),'k');
% legend([p1],'Average GCaMP');
% savefig([section_type ' ( ' num2str(cutoff_df) '% cutoff) ' num2str(timeExtension) 's no alignment std avg MAY2020 ZSCORE all syll.fig']);

% figure;
% x1 = all_bout_lengths(all_bout_lengths~=0);
% x2 = transpose(init_bout_lengths(init_bout_lengths~=0));
% x3 = first3_bout_lengths(first3_bout_lengths~=0);
% 
% x = [x1; x2; x3];
% g = [ones(size(x1)); 2*ones(size(x2)); 3*ones(size(x3))];
% boxplot(x,g);
% title([section_type ' sentence lengths (seconds)']);
% set(gca,'XTickLabel',{'All sentence lengths','Initial sentence lengths', 'First 3 sentence lengths'})
% xlabel ('');
% ylabel('seconds');

% figure;
% boxplot(init_bout_lengths);
% title('Average of first sentence lengths (seconds)')
% xlabel ('');
% ylabel('seconds');
% Saving skipped ID's/ trials to file
fid = fopen([section_type ' skipped' '.txt'],'w');
str = [repmat('%s',1,size(string(skipped_IDs),2)),'\n'];
fprintf(fid, str, transpose(string(skipped_IDs)));
fclose(fid);

% saved ID's
fid = fopen([section_type ' saved' '.txt'],'w');
str = [repmat('%s',1,size(string(saved_IDs),2)),'\n'];
fprintf(fid, str, transpose(string(saved_IDs)));
fclose(fid);


%% plot diff in average of gcamp from baseline vs section type timepoints and boxplot
% figure;
% hold on;
% % avg_gcamp=nanmean(nousv_scrambledSig_norm, 1);
% noUSV_avgSig= nanmean(no_usv_sig_norm, 2);
% noUSV_avgSig = rmmissing(noUSV_avgSig);
% spontaneousUSV_avgSig= nanmean(avg_sig_norm_cleaned, 2);
% spontaneousUSV_avgSig = rmmissing(spontaneousUSV_avgSig);
% 
% [h,p,ci,stats] = ttest2(nanmean(no_usv_sig_norm, 2),nanmean(avg_sig_norm_cleaned, 2))
% 
% title(['Z-score Avg. GCaMP baseline vs. ',[section_type],' (all syll)', ' p=', num2str(p)]);
% plot([nanmean(noUSV_avgSig), nanmean(spontaneousUSV_avgSig)],'--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
% ylabel('Z-score normalized dF/F'); % each signal was normalized using zscore before processing
% 
% % ylim([-0.05 0.25])
% savefig(['Z-score avg. gcamp nousv baseline vs ', section_type ,' (all syll).fig']);
% 
% save(['zscore gcamp_at_usv_',section_type ,'.mat'], 'avg_sig_norm_cleaned')
% 
% figure;
% G = transpose([transpose(ones(size(noUSV_avgSig)))  transpose(2*ones(size(spontaneousUSV_avgSig)))]);
% X = transpose([transpose(noUSV_avgSig), transpose(spontaneousUSV_avgSig)]);
% boxplot(X,G,'notch','on','colors',[0 0 0],'symbol','','labels',{'No-USV',section_type});
% ylabel('Z-score normalized dF/F');
% title(['Boxplot of Z-score GCaMP at Random Baseline timepoints vs. ',[section_type],' (all syll)', ' p=', num2str(p)]);
% ylim([-3, 4]);
% 
% meanDiff =  nanmean(nanmean(avg_sig_norm_cleaned, 2)) - nanmean(nanmean(no_usv_sig_norm, 2));
% dim = [.35 .5 .3 .3];
% annotation('textbox',dim,'String',['Mean Difference: ' num2str(meanDiff)],'FitBoxToText','on');
% 
% savefig(['BOXPLOT Z-score avg. gcamp baseline vs ', section_type ,' (all syll).fig']);




%% plot diff in average of gcamp from no-usv SAME SECTION vs SAME SECTION syll timepoints and boxplot
% figure;
% hold on;
% % avg_gcamp=nanmean(nousv_scrambledSig_norm, 1);
% noUSV_avgSig= nanmean(no_usv_Same_Section, 2);
% noUSV_avgSig = rmmissing(noUSV_avgSig);
% spontaneousUSV_avgSig= nanmean(avg_sig_norm_cleaned, 2);
% spontaneousUSV_avgSig = rmmissing(spontaneousUSV_avgSig);
% 
% [h,p,ci,stats] = ttest2(nanmean(no_usv_Same_Section, 2),nanmean(avg_sig_norm_cleaned, 2))
% 
% title(['Z-score Avg. GCaMP No-USV Same Section vs. ',[section_type],' (all syll)', ' p=', num2str(p)]);
% plot([nanmean(noUSV_avgSig), nanmean(spontaneousUSV_avgSig)],'--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])
% ylabel('Z-score normalized dF/F'); % each signal was normalized using zscore before processing
% 
% % ylim([-0.05 0.25])
% xlim([0 3]);
% savefig(['Z-score avg. gcamp no-usvs same section vs ', section_type ,' (all syll).fig']);
% 
% % save(['zscore gcamp_at_usv_',section_type ,'.mat'], 'avg_sig_norm_cleaned')
% 
% figure;
% G = transpose([transpose(ones(size(noUSV_avgSig)))  transpose(2*ones(size(spontaneousUSV_avgSig)))]);
% X = transpose([transpose(noUSV_avgSig), transpose(spontaneousUSV_avgSig)]);
% boxplot(X,G,'notch','off','colors',[0 0 0],'symbol','','labels',{'No-USV',section_type});
% ylabel('Z-score normalized dF/F');
% title(['Boxplot of Z-score GCaMP at Random Same Section timepoints vs. ',[section_type],' (all syll)', ' p=', num2str(p)]);
% ylim([-3, 4]);
% 
% meanDiff =  nanmean(nanmean(avg_sig_norm_cleaned, 2)) - nanmean(nanmean(no_usv_Same_Section, 2));
% dim = [.35 .5 .3 .3];
% annotation('textbox',dim,'String',['Mean Difference: ' num2str(meanDiff)],'FitBoxToText','on');
% 
% savefig(['BOXPLOT Z-score avg. gcamp same section vs ', section_type ,' (all syll).fig']);
% 
