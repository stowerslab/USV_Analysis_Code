load ('BG9-T1-FU_000_processed.mat');
GGGcamp = sig_norm(1,:);


idx_end_sentences = []; % record the end of sentences
idx_init_sentences = []; % record the begining of sentences
manipulate_timepoint = timepointUSVs;
for i=1:length(interUSVinterval)
    if interUSVinterval(i)>=1.5 && interUSVinterval(i+1) >= 1.5 % set the gap between sentences is 1.5s
        manipulate_timepoint(i) = 0;
        manipulate_timepoint(i+1) = 0;
    elseif interUSVinterval(i)>=1.5
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

% for i=1:2:length(timepointUSVs)-1
%     if timepointUSVs(i+1)-timepointUSVs(i) >= 2 % set the gap btw sentences is 2 sec
%         
%         idx_end_sentences = [idx_end_sentences, i]; % store the end time idx into the pauses array
%         idx_init_sentences = [idx_init_sentences, i+1]; % store the begin time idx into the pauses array
%         
%     end
% end
% idx_end_sentences = [idx_end_sentences,length(timepointUSVs)];
% 
% % remove bursts of noise
% %for i=1:(length(idx_init_sentences)-1)
% for i=1:length(idx_init_sentence)
%     %disp(i); 
%     if (idx_end_sentences(i+1)-idx_init_sentences(i))==1
%         idx_end_sentences(i+1)=[];
%         idx_init_sentences(i)=[];
%     end  
% end

% cut usvPowerPerSampleSmooth into sentence pieces
% cut Gcamp vector into sentence pieces
% init_sen = 1;
ave_corr = 0;
for k=1:length(idx_init_sentences)
    end_sen = idx_end_sentences(k+1); % extract the idx of last peaks in a sentence
    %disp(['end_sen: ', end_sen]);
    init_sen = idx_init_sentences(k);
    %disp(['init_sen: ', init_sen]);
    onset = locUsvs(init_sen)+20;
    offset = locUsvs(end_sen)+20;
    sen_usv = usvPowerPerSampleSmooth(onset:offset);
    %save(sen_file_usv,'sen_usv','-append');
    
    
    % calculate the time corresponding to the sentences. Using normal
    % timescle
    sen_duration_usv = timescaleUSV(onset:offset);
    %init_time = sen_duration_usv(1);
    %end_time = sen_duration_usv(length(sen_duration_usv)-1); % I am not sure if I need to save time duration for that
    init_time = timescaleUSV(onset);
    end_time = timescaleUSV(offset);
    
    % find corresponding indices from Gcamp timescale
    %init_time_Gcamp = find(timescaleGcamp == init_time*samplerateRatio);
    init_time_Gcamp = round(onset/samplerateRatio,0);
    %end_time_Gcamp = find(timescaleGcamp == end_time*samplerateRatio);
    end_time_Gcamp = round(offset/samplerateRatio,0);
    sen_duration_Gcamp = timescaleGcamp(init_time_Gcamp:end_time_Gcamp);
    
    % use such indices to find corresponding Gcamp signals
    sen_Gcamp = GGGcamp(init_time_Gcamp:end_time_Gcamp);
    %save(sen_file_Gcamp,'sen_Gcamp', '-append');
   
    init_sen = idx_init_sentences(k); % update the idx of the next peaks
    
    % downsample sentences to calculate correlation
    idx = 1:length(sen_usv);                                 % Index
    idxq = linspace(min(idx), max(idx), length(sen_Gcamp));    % Interpolation Vector
    usv_down = interp1(idx, sen_usv, idxq, 'linear');       % Downsampled Vector
    cor = corrcoef(usv_down,sen_Gcamp);
    cor = cor(2,1);
    ave_corr = ave_corr+cor;
    figure;
%     disp(['length of x: ',length(timescaleUSV(locUsvs(init_sen):locUsvs(end_sen)+1))]);
%     disp(['length of y: ',length(end_sen)]);
    plot(timescaleUSV(onset:offset), sen_usv,'r');
    plot(sen_duration_usv, sen_usv,'r');
    hold on;
    plot(sen_duration_Gcamp, sen_Gcamp, 'g');
end
ave_corr = ave_corr/length(idx_init_sentences);
