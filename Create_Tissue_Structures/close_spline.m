function coord = close_spline(coord)
    % Remove duplicate points and close a spline (first = last point)

    % Transpose the matrix to get N x 2 where each row is a point
    coord_t = coord.';
    % Use 'unique' function to find the unique rows, 'rows' option is for rows,
    % 'stable' option preserves the order of first occurrences
    [~, ia, ~] = unique(coord_t, 'rows', 'stable');
    % Extract the unique coordinates using the index array `ia`
    unique_coord_t = coord_t(ia, :);
    % Transpose back to get 2 x M matrix
    coord = unique_coord_t.';

    coord(:,end) = coord(:,1); 

end