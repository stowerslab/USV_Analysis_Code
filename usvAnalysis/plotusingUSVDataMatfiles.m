wavPath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\BNST ChR2 USV Master dataset\10Hz_merged\Male data only\';
files = dir('*.mat');

totalMice = size(files, 1);
usvPowerPerSampleSmoothMat = zeros(1,22463); 
for j = 1:totalMice
 load( files(j).name);
    usvPowerPerSampleSmoothMat = cat(1, usvPowerPerSampleSmoothMat,...
         usvPowerPerSampleSmooth); % PETE- originally concat usvPowerPerSampleSmooth

    fprintf(1, '\n');               % put a line break between each sample
   
end
usvPowerPerSampleSmoothMat = usvPowerPerSampleSmoothMat(2:end,:); 

% take the average along the columns
usvPowerPerSampleSmoothAvg = mean(usvPowerPerSampleSmoothMat);


figure;
hold on;

% plot the average of the smoothed out signals
boundedline(T(1:validSamples), mean(usvPowerPerSampleSmoothMat(:, :, 1)), std(usvPowerPerSampleSmoothAvg(:, :, 1)));

% plot(T(1:validSamples), usvPowerPerSampleSmoothAvg, 'Color', [ 0.5843 0.8157 0.9882]); % r: - PETE
title(''); % note: next one overlaps it anyway
ylabel({'USV power';'(total dB in 40-90kHz band)'});
xlabel('time (seconds)')
