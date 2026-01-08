
clc; clearvars; close all hidden;

%% Set paths and read image
impath = '/home/aferrara/Desktop/abaqus-compaction-model/Create_Tissue_Structures'; % path to image folder
imname = 'testspruce_EW3'; % image name
img = imread(fullfile(impath,strcat(imname,'.bmp'))); % read image

%% Get regions properties and filter out zero-area fibers
% Set regions properties
L = bwlabeln(img,4); % labels for connected components using 4-connectivity
stats = regionprops(L,'Area','BoundingBox','Image','Centroid','Perimeter','Conveximage','ConvexHull'); % get regions properties

out_idx = [];
for i = 1:length(stats) % for each region
    iimg = imcomplement(stats(i).Image);
    Limg = bwlabeln(iimg);
    centpos = (round(size(Limg)/2));
    lid = Limg(centpos(1),centpos(2));
    iimg(Limg~=lid) = 0;
    iimg_stats = regionprops(iimg,'Area','ConvexHull');
    try stats(i).lumenArea = iimg_stats.Area;
        stats(i).lumenConvexHull = iimg_stats.ConvexHull;
    catch; out_idx(end+1,1) = i;
    end
    
end
stats(out_idx) = []; % remove zero-area fibers

%% Fit splines
figure;
f = waitbar(0/length(stats),'Please wait...','Name','Progress'); % start wait bar
out_idx = [];
for i=1:length(stats) % for each region (=fiber)
    stats(i).skel = bwskel(stats(i).Image); % get skeleton line of object
    skelstats = regionprops(stats(i).skel,'Perimeter','Area'); % get perimeter of skeleton
    stats(i).wthick = stats(i).Area/skelstats(1).Perimeter; % mean wall thickness = area/skeleton length
    % Calsulate outer spline
    stats(i).outerpts = resample_spline(stats(i).ConvexHull);
    plot(stats(i).outerpts(1,:),stats(i).outerpts(2,:),'-b'); hold on 
    % Calculate inner spline
    innpts = resample_spline(stats(i).lumenConvexHull);
    stats(i).innerpts = innpts+[stats(i).BoundingBox(1);stats(i).BoundingBox(2)].*ones(size(innpts)); % store innerpoints
    plot(stats(i).innerpts(1,:),stats(i).innerpts(2,:),'-r'); hold on    
    % Check for intersections
    [xi, yi] = polyxpoly(stats(i).innerpts(1,:),stats(i).innerpts(2,:),...
        stats(i).outerpts(1,:),stats(i).outerpts(2,:)); % intersection coordinates
    if ~isempty(xi)
        out_idx(end+1,1) =i;
    end % store twisted fiber position
   
    waitbar(i/length(stats),f,sprintf('Please wait... %d/%d',i,length(stats)));
end
close(f)
close(gcf);
stats(out_idx) = []; % remove twisted fibers
stats(106:end) = []; % remove fibers for EW3 (right side)

% Plot fibers
figure;
for i=1:length(stats)
    plot(stats(i).outerpts(1,:),stats(i).outerpts(2,:),'-b'); hold on
    plot(stats(i).innerpts(1,:),stats(i).innerpts(2,:),'-r'); hold on
    text(stats(i).Centroid(1), stats(i).Centroid(2), num2str(i), 'FontSize', 8, 'Color', 'b');
end

%% Fix fibers intersections
for i = 1:length(stats)
    for k = 1:i-1 % for each previous region
        [xi, yi] = polyxpoly(stats(i).outerpts(1,:),stats(i).outerpts(2,:),...
            stats(k).outerpts(1,:),stats(k).outerpts(2,:));
        if ~isempty(xi)
            plot(xi(:), yi(:), 'g--'); scatter(xi, yi, 'g*'); hold on
            % Replace intersection with linear points
            fiber = [i k];
            for rep=1:2
                f = fiber(rep);
                [~,rowidx1] = min(vecnorm((stats(f).outerpts'-[xi(1) yi(1)]).')); % point closest to int1
                [~,rowidx2] = min(vecnorm((stats(f).outerpts'-[xi(2) yi(2)]).')); % point closest to int2
                [~,rowidxm] = min(vecnorm((stats(f).outerpts'-[(xi(1)+xi(2))/2 (yi(1)+yi(2))/2]).')); % point closest to mean int
                if rowidx1~=rowidx2
                    if ~ismember(rowidxm, [rowidx1,rowidx2])
                        % Create array of indexes between intersection points
                        if rowidxm<max(rowidx1,rowidx2) && rowidxm>min(rowidx1,rowidx2) % if indexes contiuous between intersection points
                            if rowidx1<rowidx2; idx = rowidx1:rowidx2;
                            else; idx = rowidx1:-1:rowidx2; end
                        else % if indexes non continuous between intersection points
                            if rowidx1<rowidx2; idx = [rowidx1:-1:1, length(stats(f).outerpts):-1:rowidx2];
                            else; idx = [rowidx1:length(stats(f).outerpts), 1:rowidx2]; end
                        end
                    else
                        idx = [rowidx1, rowidx2];
                    end
                    % Replace y with linear values
                    stats(f).outerpts(1,idx(:)) = linspace(xi(1),xi(2),length(idx)); % replace x with equally spaced value between intersection points
                    m = (yi(2)-yi(1))/(xi(2)-xi(1)); % slope
                    b = yi(1)-m*xi(1); % intercept
                    
                    stats(f).outerpts(2,idx(:)) = m*stats(f).outerpts(1,idx(:))+b; % replace y
                    % Adjust fibers for EW3
                    if i==94 && k==77 && rep==1
                        stats(f).outerpts(:,[idx(1),idx(end)]) = [];             
                    end
                    plot(stats(f).outerpts(1,idx(:)), stats(f).outerpts(2,idx(:)), '--c'); hold on
                end
            end
            % Check if there are intersection with inner spline
            [xii,yii] = polyxpoly(stats(i).outerpts(1,:),stats(i).outerpts(2,:),...
                stats(k).innerpts(1,:),stats(k).innerpts(2,:));
            if ~isempty(xii)
                scatter(xii,yii,'c*'); hold on
            end       
        end
    end
end

%% Scale fibers to real dimension
scalepxl = 1/9;
% Scale interface structure
figure;
downscale = 0.995;
layer_thick = [0.175, 0.175, 0.125, 0.125, 0.035];
for i = 1:length(stats)

    if i == 10 || i == 101
        downscale = 0.97;
    end

    % Scale to real dimensions (micrometer)
    stats(i).Centroid = stats(i).Centroid * scalepxl;
    stats(i).outerpts = stats(i).outerpts * scalepxl;
    stats(i).wthick = stats(i).wthick * scalepxl;
    S2thick = stats(i).wthick*downscale - sum(layer_thick);
    stats(i).lthick = [layer_thick(1:end-1), S2thick, layer_thick(end)];
    % Scale outer spline to get inner spline
    stats(i).innerpts = close_spline(scale_spline(stats(i).outerpts, stats(i).Centroid, stats(i).wthick));
    % Down-scale outer spline    
    stats(i).outerpts = close_spline(scale_spline(stats(i).outerpts, stats(i).Centroid, stats(i).wthick *(1-downscale)));    
    stats(i).wthick_scaled = stats(i).wthick *downscale;
    stats(i).wthick_diff = stats(i).wthick - stats(i).wthick_scaled;

    % Compute Horizontal & Vertical Widths
    x_coords = stats(i).outerpts(1, :); % X-coordinates of outer perimeter
    y_coords = stats(i).outerpts(2, :); % Y-coordinates of outer perimeter
    % Find horizontal width: Distance between leftmost and rightmost points
    left_x = min(x_coords);
    right_x = max(x_coords);
    horiz_widths(i) = abs(right_x - left_x);
    % Find vertical width: Distance between topmost and bottommost points
    top_y = min(y_coords);  % Lowest Y (image top)
    bottom_y = max(y_coords);  % Highest Y (image bottom)
    vert_widths(i) = abs(bottom_y - top_y);

    % Plot fibers
    plot(stats(i).outerpts(1,:),stats(i).outerpts(2,:),'-b'); hold on
    plot(stats(i).innerpts(1,:),stats(i).innerpts(2,:),'-r'); hold on
    text(stats(i).Centroid(1), stats(i).Centroid(2), num2str(i), 'FontSize', 12, 'Color', 'b'); hold on
end

% Compute mean horizontal and vertical width
mean_horiz_width = mean(horiz_widths);
mean_vert_width = mean(vert_widths);

% Display results
disp(['Mean cell wall thickness: ', num2str(mean([stats.wthick]), '%.3f'), ' μm']);
disp(['Mean S2 thickness: ', num2str(mean(arrayfun(@(s) s.lthick(5), stats)), '%.3f'), ' μm']);
disp(['Mean tangential width: ', num2str(mean_horiz_width, '%.3f'), ' μm']);
disp(['Mean radial width: ', num2str(mean_vert_width, '%.3f'), ' μm']);



%% Create ligning mask
% Store outer points of all splines
all_outerpts = [];
for i=1:length(stats)
    all_outerpts = [all_outerpts; stats(i).outerpts'];
end
all_outerpts = unique(all_outerpts, 'rows', 'stable');
bound = boundary(all_outerpts(:,1), all_outerpts(:,2), 0.9);
maskpts = all_outerpts(bound,:)';
% Scale mask to include ML
C = calculate_centroid(maskpts);
maskpts = close_spline(scale_spline(maskpts, C, -3*layer_thick(1)));
% Uncomment only to find the following idx by picking points on the plot
% (read coordinates and put values below)
%figure; plot3(maskpts(1,:), maskpts(2,:), [1:length(maskpts)], 'k'); view([0, 90]);
% Flatten top and bottom matrix
b = [max(maskpts(2,:)), min(maskpts(2,:))]; % intercept
idx_left_top = 2074;
idx_right_top = 1279;
idx_left_bottom = 6;
idx_right_bottom = 1110;
ref_left_top = maskpts(2,idx_left_top);
ref_right_top = maskpts(2,idx_right_top);
idx_arr = {idx_right_top:idx_left_top, idx_left_bottom:idx_right_bottom};
for i=1:2
    idx = idx_arr{i};
    x = [maskpts(1,idx(1)), maskpts(1,idx(end))];
    % Replace y with linear values
    maskpts(1,idx(:)) = linspace(x(1),x(2),length(idx)); % replace x with equally spaced value between intersection points
    m = 0; % slope
    maskpts(2,idx(:)) = m*maskpts(1,idx(:))+b(i); % replace y
end
% Add vertical points towards top
n = 10;
% Top left
y_step = (max(maskpts(2,:))-ref_left_top)/n;
x_coords = maskpts(1,idx_left_top) * ones(1, n);
y_coords = (ref_left_top + (n-1) * y_step) : -y_step : ref_left_top;
add_coord = [x_coords; y_coords];
before_insertion = maskpts(:, 1:idx_left_top);
after_insertion = maskpts(:, idx_left_top+1:end);
maskpts = [before_insertion, add_coord, after_insertion];
% Top right
y_step = (max(maskpts(2,:))-ref_right_top)/n;
x_coords = maskpts(1,idx_right_top) * ones(1, n);
y_coords = ref_right_top : y_step : (ref_right_top + (n-1) * y_step);
add_coord = [x_coords; y_coords];
before_insertion = maskpts(:, 1:idx_right_top-1);
after_insertion = maskpts(:, idx_right_top:end);
maskpts = [before_insertion, add_coord, after_insertion];

% Plot mask
plot(maskpts(1,:), maskpts(2,:), 'k');

%% Save data for Abaqus
mstruct = struct();
mstruct.mask = maskpts;
for i=1:length(stats)
    % Store data
    mstruct.fibers(i).outerpts = double(stats(i).outerpts);
    mstruct.fibers(i).innerpts = double(stats(i).innerpts);  
    mstruct.fibers(i).Centroid = double(stats(i).Centroid);
    mstruct.fibers(i).wthick = double(stats(i).wthick); 
    mstruct.fibers(i).lthicks = double(stats(i).lthick);
end

% Save data struct
savefile = fullfile(impath, strcat(imname,'_mstruct.mat')); 
save(savefile,'mstruct');
