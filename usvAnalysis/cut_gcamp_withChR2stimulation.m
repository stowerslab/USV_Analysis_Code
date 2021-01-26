%% Cutting Gcamp into sections to match wav sections
clear all; close all;
load ('CGC4_T4_000_000_002_processed');
micename = 'CGC4_T4_';
ROW_IN_EXCEL = 5; % need to specify which row in the excel timing sheet goes
% this gcamp data

GGGcamp = sig_norm(1,:); % load one side of gcamp signal

Fs = 20; % gcamp samples at 20hz
fontSz = 25;
stimLength = 5; %in seconds
preSecToPlot = 10;
postSecToPlot = 20;

% %For 10s stimulation long
% stimLength = 5; %in seconds
% preSecToPlot = 10;
% postSecToPlot = 20;
FilePath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\GCaMP Data\2020-6 LPOA ChR2 PAG FP\2020-7-2 ChR2 PAG FP_T4 ChR2\';
wavPath = '\\172.29.164.29\Jingyi Chen\Dropbox\Jingyi Data backup\Behavior data\2020-6 chr2 PAG-FP\2020-7-3 ChR2 T5\';
xlFilePath = [FilePath 'Timesheet_T4'];
timestampData = xlsread(xlFilePath); % read in timing to split gcamp into sections

gcampStart = timestampData(ROW_IN_EXCEL,2); % where gcamp starts in the video

all_timestamps = (timestampData(ROW_IN_EXCEL,3:end)) - gcampStart; % grabs row in excel sheet with video timing

% and subtract the time where gcamp isnt being recorded
for i = [11 12 13]
   cutstart = all_timestamps(i)-preSecToPlot;
   cutend = all_timestamps(i)+postSecToPlot;
   section = round([cutstart cutend] * Fs); 
  GGGcampTrim = GGGcamp(section(1):section(2));
  matName= [FilePath, micename, 'Stimulation', int2str(i), '.mat'];
save(matName,'GGGcampTrim');
end
    