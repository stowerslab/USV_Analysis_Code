%%plot GCaMP with ChR2 stimulation
clear all; close all;
%plot boundline figures
gcampFrameRate=20;
timeExtension = 10;
GCaMP_Normed = []
files = dir('*.mat');
totalMice = size(files, 1); 

for j = 1:totalMice                     
load (files(j).name);
GCaMP_Normed (j,:) = GGGcampTrim *100;
end

figure;
hold on;
ylabel ('dF/F (%)');
xlabel ('time (s)');
% p1 = plot(-timeExtension:1/gcampFrameRate:timeExtension,mean(avg_sig_norm_cleaned, 1),'g');
p1 = boundedline(-10:1/gcampFrameRate:20,nanmean(GCaMP_Normed, 1), nanstd(GCaMP_Normed)/sqrt(size(GCaMP_Normed, 1)),'g'); % standard error of the mean is usually estimated as the sample standard deviation divided by the square root of the sample size
savefig('with 50Hz ChR2 stimulation.fig');

% ZSCORE SCALE
figure;
hold on;
ylabel ('Z-score');
xlabel ('time (s)');
% p1 = plot(-timeExtension:1/gcampFrameRate:timeExtension,mean(avg_sig_norm_cleaned, 1),'g');
p1 = boundedline(-10:1/gcampFrameRate:20,zscore(nanmean(GCaMP_Normed, 1)), nanstd(GCaMP_Normed)/sqrt(size(GCaMP_Normed, 1)),'g'); % standard error of the mean is usually estimated as the sample standard deviation divided by the square root of the sample size
% p2 = plot(-timeExtension:1/gcampFrameRate:timeExtension,zscore(nanmean(avg_scrambledSig_norm, 1)),'k');
legend([p1],'Average GCaMP');
savefig('Z score with 50Hz ChR2 stimulation.fig');