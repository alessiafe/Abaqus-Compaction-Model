clc; clear; close all;

impath = '/home/aferrara/Desktop/abaqus-compaction-model/Create_Tissue_Structures';
imname = 'testspruce_LW3';
img = imread(fullfile(impath,strcat(imname,'.bmp')));

figure;
imshow(img);
title('Draw lines by clicking two points per line. Close the figure to stop.');
hold on;

angles_rad = []; % Store angles in radians
angles_deg = []; % Store angles in degrees
horizontal_lengths = []; % Store horizontal components

disp('Start drawing lines. Select two points per line.');
disp('Close the figure to stop.');

while true
    try
        % Get two points from user
        [x, y] = ginput(2);
        
        % If user closes the figure, exit the loop
        if isempty(x) || length(x) < 2
            break;
        end

        % Draw the line
        plot(x, y, 'ro-', 'MarkerSize', 10, 'LineWidth', 2);

        % Compute the angle (in radians)
        dx = x(2) - x(1); % Horizontal difference
        dy = y(1) - y(2); % Vertical difference (image coordinates are flipped)
        angle_rad = atan2(dy, dx); % Angle in radians

        % Ensure all angles are positive (convert -π to π, etc.)
        if angle_rad < 0
            angle_rad = angle_rad + pi;
        end

        % Convert to degrees
        angle_deg = rad2deg(angle_rad);

        % Compute the horizontal component
        horizontal_length = abs(dx);

        % Store values
        angles_rad = [angles_rad, angle_rad];
        angles_deg = [angles_deg, angle_deg];
        horizontal_lengths = [horizontal_lengths, horizontal_length];

        % Display current angle
        disp(['Line ', num2str(length(angles_rad)), ': Angle = ', num2str(angle_rad, '%.3f'), ' rad (', num2str(angle_deg, '%.3f'), '°)']);

    catch
        break; % If figure is closed, exit loop
    end
end

% Compute average angle
if ~isempty(angles_rad)
    avg_angle_rad = mean(angles_rad);
    avg_angle_deg = mean(angles_deg);
    
    disp('----------------------------------');
    disp(['All angles (radians): ', num2str(angles_rad, '%.3f ')]);
    disp(['All angles (degrees): ', num2str(angles_deg, '%.3f ')]);
    disp(['Average angle: ', num2str(avg_angle_rad, '%.3f'), ' rad (', num2str(avg_angle_deg, '%.3f'), '°)']);
end
