function pts_new = scale_spline(pts, centroid, thick)
    % Scale a spline by a given thickness (outward -, inward +), 
    % keeping the same centroid, so that the new spline lies at a constant
    % distance from the original.
    pts_new = []; % initialize array
    % Translate outer spline
    trans_x = pts(1,:) - centroid(1);
    trans_y = pts(2,:) - centroid(2);
    % Scale spline
    scalef = 1 - (thick / sqrt(mean(trans_x.^2 + trans_y.^2)));
    pts_new(1,:) = centroid(1) + scalef * trans_x;
    pts_new(2,:) = centroid(2) + scalef * trans_y;

end