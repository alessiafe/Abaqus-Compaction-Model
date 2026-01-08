function pts = resample_spline(pts)
    % Re-sample spline to increase points density for better distribution
    pts = cscvn(pts');
    tout = ppval(pts,linspace(0,pts.breaks(end),2000)); % resampling for better point distribution
    dtout = diff(tout,1,2);
    db = pts.breaks(end)/500;
    fout = tout(:,1); dist = 0;
    for zi = 1:length(dtout)-1
        dist = dist+norm(dtout(:,zi));
        if dist >= db; fout = [fout tout(:,zi)]; dist=0; end
    end
    fout=[fout tout(:,end)];
    pts = fout; % store outerpoints

end