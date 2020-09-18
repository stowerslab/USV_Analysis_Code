clear all; close all;

matName = 'C:\data\Jason\microscope\2016_01_CrhTdtwithEsrImmuno\male\testBar.mat';
load(matName);
[thresh, IbwNissl, nisslCC] = sectionSegmentCountRoi(IsubNissl);

nisslCC.NumObjects

