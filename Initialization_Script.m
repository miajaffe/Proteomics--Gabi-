% Initialization Script
% Run me!
clear all
close all
clc
filename = 'Longitudinal_RawCounts_ForClass.xlsx';
[OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename);