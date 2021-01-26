%JC for combine data from raw .mat files and fit timepointUSVs distribution for
%number of USV plot across time
clear all; close all;
wavPath = '\\172.29.164.29\jingyi chen\Dropbox\Jingyi Data backup\BNST ChR2 USV Master dataset\PAG VGlu USVs\T1T2\5Hz\';
%Note have to put allstats.mat into another foder
files = dir('*.mat');
totalMice = size(files, 1);
nameList = {};
a = 0;
%first combine all the data into one big matrix
for j = 1:totalMice
   nameList = [nameList files(j).name];
   load (files(j).name);
   start=a;
   finish = start + length (timepointUSVs);
   Timepointcombined(start+1:finish,1) = timepointUSVs';
   a=finish;
end
%now plot distrubution and try to fit with kernel
[f,c] = hist(Timepointcombined,100); %bin size: 100, can change based on data size
% normalization to get density
f = f/trapz(c,f);
% kernel density
pts = linspace(0, 20, 100);
[fk,xk] = ksdensity(Timepointcombined, pts, 'Support', [0, 20]);
figure;
bar(c,f,'w'); % plot the distrubution in bars
hold on
plot(xk,fk, 'r', 'LineWidth', 2); % plot the kernel fit in red
% ylabel({'Number of USV'}); %frequency of the number count, what is the unit?
xlabel('time (seconds)');
xlim([0 20]); % plot x axis 
set(gca, 'FontSize', 20);     % font size of 20
%plot the stimulus box
xLims = get(gca, 'XLim');
axis([xLims(1) xLims(2) -0.15 1]);    % -0.05 to go slightly below lowest call, 3 should be above highest call for all VGlu T1T2 ChR2 and female ChR2
% axis([xLims(1) xLims(2) -0.15 1.2]);     % -0.05 to go slightly below lowest call, 2 should be above highest call for all other ChR2
% axis([xLims(1) xLims(2) -0.15 3.5]);  % -0.05 to go slightly below lowest call, 2.5 should be above highest call for all VGlu T3 ChR2
x = [5 10 10 5];
yLims = get(gca, 'YLim');
y = [yLims(1) yLims(1) yLims(2) yLims(2)];
% patch(x, y, [1 0.3 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % orange box
patch(x, y, [0 0.4 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % blue BOX FOR STIMULATION PERIOD
hold off;

%%Use hisfit to fit the distribution
figure;
hold on;
h=histfit(Timepointcombined,100,'kernel');% use 100 bins to make the plot looks better
h(1).FaceColor = [.8 .8 1];
h(2).Color = [0.9100 0.4100 0.1700];
ylabel({'Number of USV'});
xlabel('time (seconds)');
xlim([0 20]); % plot x axis 
set(gca, 'FontSize', 20);     % font size of 20
hold off;


