function C = calculate_centroid(pts)
    % Calculate centroid coordinates of a spline
    
    % Calculate the area of the polygon using the shoelace formula
    A = 0.5*abs(sum(pts(1,1:end-1).*pts(2,2:end)-pts(1,2:end).*pts(2,1:end-1)));
    % Calculate the centroid coordinates
    C(1) = (1/(6*A))*sum((pts(1,1:end-1)+pts(1,2:end)).*(pts(1,1:end-1).*pts(2,2:end)-pts(1,2:end).*pts(2,1:end-1)));
    C(2) = (1/(6*A))*sum((pts(2,1:end-1)+pts(2,2:end)).*(pts(1,1:end-1).*pts(2,2:end)-pts(1,2:end).*pts(2,1:end-1)));

end