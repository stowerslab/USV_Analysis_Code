% top level script to count all VMH sections and plot results
clear all; close all;

% %TMP-L virus isoflurane
% filepathRoot = 'C:\data\jason\peterVirusData\vmhIso';  %root for all mice
% mouseNums = {'v23'};  %ex. 2nd entry is mouseNums{1,2}
    
totalMice = size(mouseNums, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now uncomment below to save JPEGs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    
    matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'] %#ok<NOPTS>
    load(matName);
    %RED:
    hf = figure; imagesc(Ir); %use 'imagesc' linearly scale to full bit resolution 
    colormap(gray); axis off;
    set(hf, 'Position', get(0,'Screensize')); % maximize figure
    % clear the figure background before saving:
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')
    saveas(hf,[filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_vlTdt.jpg'], 'jpg');% 
    close(hf); %prevent large number figures
    %GREEN:
    hf = figure; imagesc(Ig); %use 'imagesc' linearly scale to full bit resolution 
    colormap(gray); axis off;
    set(hf, 'Position', get(0,'Screensize')); % maximize figure
    % clear the figure background before saving:
    whitebg('black')
    set(gcf,'Color',[0 0 0])
    set(gcf,'InvertHardcopy','off')
    saveas(hf,[filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_vlGfp.jpg'], 'jpg');% 
    close(hf); %prevent large number figures 
        
 
end




