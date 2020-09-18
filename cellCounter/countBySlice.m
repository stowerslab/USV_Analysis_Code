clear;

% put in the correct channel#CellCoord variable at line 14

countsBySlice = [];

sectionFilePath = 'E:\Sourish\FOR PETE'; % select directory with the .mat files
sectionDirectory = dir(fullfile(sectionFilePath,'*.mat'));

for sectionData = sectionDirectory'
    load(fullfile(sectionData.folder, sectionData.name)); 
    sectionData.name % need to make sure files are in order numerically,
    %can manually change digits in filenames to have a leading zero
    countsBySlice = [countsBySlice, length(channel3CellCoord)]; % PUT IN THE CORRECT CHANNEL
end
rawCounts = countsBySlice;
countsBySlice = countsBySlice/max(countsBySlice);
plot(linspace(1, length(countsBySlice), length(countsBySlice)), countsBySlice, '-o')
xlabel('slice #')
ylabel('normalized cell count')
title('Infection through brain slices')