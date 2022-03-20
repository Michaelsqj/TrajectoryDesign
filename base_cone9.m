function [time, g, k] = base_cone9(mat, fov, Ts, gmax, smax, readout_time, cone_angle)
    
    %----------------------------
    % first design the perfusion region
    %   then move to the kmax
    
    %-------------------------
    %   central k-space
    res         = fov / mat; %mm
    kmax        = 10 / (2 * res);  % [1/cm]

    [time1, gT1, kT1] = base_cone3(ceil(mat/3), fov, Ts, gmax, smax, readout_time/3, cone_angle);

    [time2, gT2, kT2] = base_cone3(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
    
    kT1endZ = kT1(end,3);
    pt2 = find(kT2(:,3)>=kT1(end,3), 1);

    dphi = atan2(kT2(pt2,2), kT2(pt2,1)) - atan2(kT1(end,2), kT1(end,1));
    R_polar = [cos(dphi) -sin(dphi) 0;
                sin(dphi)  cos(dphi) 0;
                0          0        1];

    kT1 = (R_polar * kT1')';
    kT = [kT1; kT2(pt2:end,:)];

    % time = Ts*length(k);
    % g = [0,0,0; diff(k)/Ts/4.258];
    midpt = length(kT1);
    difkT = diff(kT);
    pts = [1:(midpt-10),(midpt+10):length(difkT)];
    midps = midpt-9 : midpt+9;
    midv(:,1) = spline(pts, [difkT(1:(midpt-10),1);difkT((midpt+10):length(difkT),1)], midps);
    midv(:,2) = spline(pts, [difkT(1:(midpt-10),2);difkT((midpt+10):length(difkT),2)], midps);
    midv(:,3) = spline(pts, [difkT(1:(midpt-10),3);difkT((midpt+10):length(difkT),3)], midps);
    
    difkT = [difkT(1:(midpt-10),:); midv; difkT((midpt+10):length(difkT),:)];
    
    kT = cumsum([0,0,0; difkT]);
    
    [C, time, g, s, k, phi, sta, stb] = minTimeGradient(kT, 0, 0, 0, gmax, smax, Ts);

    [grd, krd, lenrd] = minRampDown(g(end,:),k(end,:), 1e3*smax,1e-3*Ts);
    g = [g; grd];
    k = [k; krd];
    time = time + Ts*lenrd;
    % b = 1;
    
    % kTr1end = sqrt(sum(kT1(end,1:2).^2));
    % phi1end = atan2(kT1(end,2), kT1(end,1));
    
    % kmax_circ = 3^b*kTr1end;
    % kmax_rad = sqrt( kmax^2 - kmax_circ^2 );
    
    % t = linspace(0,1,1001);
    % mid = 1 + ceil(kT1(end,3) / kmax_rad / (t(2)-t(1)));
    % t2 = t(mid:end);
    
    % fact = 0.0;
    % p = 2;
    % for step = [1 0.5 0.1 0.05 0.01 0.005]
    %     rtime = 0;

    %     while rtime < readout_time
    %         fact = fact + step;
            

            
    %         kTr2 = t2.^b * kmax_circ;

    %         phi2 = 2 * pi * fact * (t2.^p - t2(1).^p) + phi1end;
            
    %         kT2(:,1) = kTr2 .* cos(phi2);
    %         kT2(:,2) = kTr2 .* sin(phi2);
                 

    %         kT2(:, 3) = t2.^1.0 * kmax_rad;

    %         kT = [kT1; kT2];

    %         [C, time, g, s, k, phi, sta, stb] = minTimeGradient(kT, 0, 0, 0, gmax, smax, Ts);

    %         rtime = max(time);
    %     end
    %     fact = fact - step;
    %     fprintf(['\n\n A is ',num2str(fact),'\n']);
    % end

end