function [time, g, k] = base_cone8(mat, fov, Ts, gmax, smax, readout_time, cone_angle)
    % connect the johnson's + gurney's
    res         = fov / mat; %mm
    kmax        = 10 / (2 * res);  % [1/cm]
    [time2, tmpg2, tmpk2] = base_cone3(mat, fov, Ts, gmax, smax, 0.95 * readout_time, cone_angle);
    pts = find( sqrt(sum(tmpk2.^2,2)) >= (kmax/3) );
    extptsNum = 20;
    extpts = (pts(1)-extptsNum): (pts(1)-1);
    g2 = tmpg2(pts,:);
    k2 = tmpk2(pts,:);
    extg2 = tmpg2(extpts,:);
    extk2 = tmpk2(extpts,:);

    % connect k-space center to extk2
    t = linspace(0,1,1000);
    kTrmax = sqrt( extk2(1,1).^2 + extk2(1,2).^2 );
    kTzmax = extk2(1,3);

    phifin = atan2( extk2(1,2), extk2(1,1) );

    gfin = sqrt( sum(k2(1,:).^2) );
    fact = 0.0;
    p = 1;
    for step = [1 0.5 0.1 0.05 0.01 0.005]
        rtime = 0;
        while rtime < (readout_time - Ts * size(k2,1))
            fact = fact + step;

            % the azimuthal angle of the last point of k1, should be same as 
            %   starting point of k2
            phi0 = phifin - 2 * pi * fact;     
            phi = 2 * pi * fact * t.^p + phi0;

            kT(:, 1) = cos(phi) .* t .* kTrmax;
            kT(:, 2) = sin(phi) .* t .* kTrmax;
            kT(:, 3) = t * kTzmax;

            kT = [kT; extk2];

            [C, time, g1, s, k1, phi, sta, stb] = minTimeGradient(kT, 0, 0, 0.5 * gfin, gmax, smax, Ts);

            rtime = max(time);
            clear kT
        end

        fact = fact - step;
    end
    g = [g1; g2];
    k = [k1; k2];
    assert( sum((g2(1,:) - g1(end,:)).^2) < ((Ts * smax)^2) );

end