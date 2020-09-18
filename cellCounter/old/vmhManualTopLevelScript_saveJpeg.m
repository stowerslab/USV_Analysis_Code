% top level script to count all VMH sections and plot results
clear all; close all;

% %TMP-L virus isoflurane
% filepathRoot = 'C:\data\jason\peterVirusData\vmhIso';  %root for all mice
% mouseNums = {'v5','v6','v7','v8','v12','v15','v16','v17','v20',... % INPUT folder names for individual mice:
%              'v21','v23','v24','v25','v26','v28','v29','v30'};  %ex. 2nd entry is mouseNums{1,2}

% %TMP-DMSO virus
% filepathRoot = 'C:\data\jason\peterVirusData\firstAggrTmpDmso\peterCountApril';  %root for all mice
% mouseNums = {'a1','a3','a4','c1','c2','c3','c4'}; % INPUT folder names for individual mice; %ex. 2nd entry is mouseNums{1,2}
         
%TMP-DMSO Ai9
% filepathRoot = 'C:\data\jason\peterAi9data\2014_Ai9FearAggr\aggrCastrateWithUrine';  %root for all mice
filepathRoot = 'C:\data\jason\peterAi9data\2014_Ai9FearAggr\aggrHomecage'; 
% mouseNums = {'3a','3b','3c','3d','3f'}; % INPUT folder names for individual mice        
mouseNums = {'1b','1c','1d'}; % INPUT folder names for individual mice  
    
totalMice = size(mouseNums, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do all manual steps for each section first (i.e. ROIs and thresholds):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:4  %loop over sections for each mouse (must have 4 sections per mouse numbered 1-4)
        % do image loading and ROIs for each section first:
        useBlue = true;  %true if using blue
        [Ir, Ig, Ib, hVmhVlPatchMask, ellip] = vmhSetRoiEllipse(filepathRoot, mouseNum, num2str(j), useBlue); 
        %save data for each individual section:
        matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'];
        save(matName, 'Ir', 'Ig', 'Ib', 'hVmhVlPatchMask', 'ellip');
        close all; %prevent large number figures
        
        % Use the green channel to draw ROI alternatively:
        choice = questdlg('Using green to draw another ROI?','a query?', 'Yes','No','No');
        % Handle response
        switch choice
            case 'Yes'; useBlue = false;
            case 'No'; useBlue = true;
        end
        
        if ~useBlue %make an ROI using GFP channel and overwrite that using blue
            [Ir, Ig, Ib, hVmhVlPatchMask, ellip] = vmhSetRoiEllipse(filepathRoot, mouseNum, num2str(j), useBlue); 
            %save data for each individual section:
            matName = [filepathRoot, '\', mouseNum, '\', mouseNum, '_', num2str(j), '_roiImages.mat'];
            save(matName, 'Ir', 'Ig', 'Ib', 'hVmhVlPatchMask', 'ellip');
            close all; %prevent large number figures
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now uncomment below to save JPEGs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:totalMice
    mouseNum = mouseNums{1,k};
    for j = 1:4
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
 
end




