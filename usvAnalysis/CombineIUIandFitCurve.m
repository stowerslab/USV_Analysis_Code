%JC for combine data from raw .mat files and fit IUI value distribution for
%determine the sentence breaks
clear all; close all;
wavPath = '\\172.29.164.29\jingyi chen\Dropbox\Jingyi Data backup\Behavior data\WT behaviors\2018-2 WT USVs\2018-2-28 WT USV FEMALE URINE\2MINS TRIMED AUDIO\';
%Note have to put allstats.mat into another foder
files = dir('*.mat');
totalMice = size(files, 1);
nameList = {};
a = 0;
for j = 1:totalMice
   nameList = [nameList files(j).name];
   load (files(j).name);
   start=a;
   finish = start + length (interUSVinterval);
   IUIcombined(start+1:finish,1) = interUSVinterval';
   a=finish;
end
%now plot all raw IUI value on one plot raw data to check distribution
his=histogram(IUIcombined);
saveas(his, [wavPath, 'IUIRawDistribution.tif'], 'tif');
% Now fit data based on histogram check. Feels like Poission for raw IUIs
pd=fitdist (IUIcombined, 'Poisson');
h=chi2gof(IUIcombined,'CDF',pd);

matName = [wavPath,'IUIcombined.mat']; %save data for each individual mouse
save(matName,'IUIcombined');
% pd=fitdist (IUIcombined, 'Normal');
% [h,p]=chi2gof(IUIcombined,'CDF',pd);