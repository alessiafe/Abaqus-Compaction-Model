clc; clearvars; close all hidden;

%% Set paths
impath = '/home/aferrara/Desktop/abaqus-compaction-model/Create_Tissue_Structures';
imname = 'testspruce_EW3';

%% Load file
savefile = fullfile(impath, strcat(imname, '_mstruct.mat'));
load(savefile, 'mstruct');

%% Extract data and calculate
% Extract data from the structure
stats = mstruct.fibers; % Assign to stats for easier handling
% Compute mean values
mean_wthick = mean([stats.wthick]);  % Mean cell wall thickness
mean_S2thick = mean(arrayfun(@(s) s.lthicks(5), stats));  % Mean S2 thickness
% Compute horizontal & vertical widths
horiz_widths = arrayfun(@(s) max(s.outerpts(1, :)) - min(s.outerpts(1, :)), stats);
vert_widths = arrayfun(@(s) max(s.outerpts(2, :)) - min(s.outerpts(2, :)), stats);
mean_horiz_width = mean(horiz_widths);
mean_vert_width = mean(vert_widths);
% Display results
disp(['Mean cell wall thickness: ', num2str(mean_wthick, '%.3f'), ' μm']);
disp(['Mean S2 thickness: ', num2str(mean_S2thick, '%.3f'), ' μm']);
disp(['Mean tangential width: ', num2str(mean_horiz_width, '%.3f'), ' μm']);
disp(['Mean radial width: ', num2str(mean_vert_width, '%.3f'), ' μm']);