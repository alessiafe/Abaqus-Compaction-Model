
clc; clearvars; close all hidden;

%% Set paths and read image
impath = '/home/aferrara/Desktop/abaqus-compaction-model/Create_Tissue_Structures';
basename = 'testspruce_';
tissname = {'EW3', 'LW3', 'TW3'};

mean_size = [];
for i = 1:length(tissname)
    imname = strcat(basename,tissname{i});
    filepath = fullfile(impath, strcat(imname,'_mstruct.mat'));
    data = load(filepath).mstruct;
    max_size = [];
    for k = 1:length(data.fibers)
        x = data.fibers(k).outerpts(1,:);  % x-coord        
        max_size(k) = max(x) - min(x);
    end
    mean_size(i) = mean(max_size);
end

tissname{end+1} = 'Average';
mean_size(end+1) = mean(mean_size);

% Create table
T = table(tissname', mean_size', 'VariableNames', {'Tissue', 'Mean size'});
disp(T);



