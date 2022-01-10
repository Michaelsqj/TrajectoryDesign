function [time, g, k] = base_cone7(mat, fov, Ts, gmax, smax, readout_time, cone_angle)
    % ---------------------------------------------------
    % code to connect johnson's and gurney's cone
    %   1. define the upper half shape, kTr
    %   2. phi = a*t, z = b * t, kTr = r^2
    
    res         = fov / mat; %mm
    kmax        = 10 / (2 * res);  % [1/cm]
    
    t = linspace(0, 1, 1000); % t is the 'alpha' in eq[2]
    turn_pt = 750;
    t1 = t(1:turn_pt);
    t2 = t(turn_pt + 1 : end);

    kT2(:, 3) = t2 * kmax * cos(cone_angle / 180 * pi);
    kTr2 = t2.^3 * kmax * sin(cone_angle / 180 * pi);

    [time1, tmpg1, tmpk1] = base_cone3(mat, fov, Ts, gmax, smax, readout_time, atan2d( kTr2(1), kT2(1,3)));
%     theta = atan2(kT2(1,3), kTr2(1));
% 
%     DCF = 1.0;
%     OS = 1;
%     LEN = readout_time / (1000*Ts);
%     NINT = 1;
%     addpath('./iGurney');
%     [tmpg1,tmpk1] = findcone(fov/mat,fov/10,LEN,theta,0.1,Ts/1000,1,smax * 1000,gmax,DCF);
%     [tmpg1,tmpk1,lenro] = gencone(fov/mat, fov/10, NINT, theta, 5e3, Ts/1000, smax * 1000, gmax, OS, DCF);
%     rmpath('./iGurney');

    k1_pts = find(tmpk1(:,3) < kT2(1,3));
    k1 = tmpk1(k1_pts,:);
    g1 = tmpg1(k1_pts,:);

    size(k1)
    len1 = size(k1, 1);
    time2 = readout_time - len1*Ts;

    phi0 = atan2( k1(end,2), k1(end,1) );
    
    phi0*180/(2*pi)

    fact = 0.0;
    p = 1;
    for step = [1 0.5 0.1 0.05 0.01 0.005]
        rtime = 0;

        while rtime < time2
            fact = fact + step;

            phi = 2 * pi * fact * (t2-t2(1)).^p + phi0;

            kT2(:, 1) = cos(phi) .* kTr2;
            kT2(:, 2) = sin(phi) .* kTr2;

%             kT(:,1) = [k1(:,1);kT2(:,1)];
%             kT(:,2) = [k1(:,2);kT2(:,2)];
%             kT(:,3) = [k1(:,3);kT2(:,3)];



            [C, time, g2, s, k2, phi, sta, stb] = minTimeGradient(kT2, 0, sqrt(sum(g1(end,:).^2)), 0, gmax, smax, Ts);

            rtime = max(time);
        end

        fact = fact - step;
    end
    fprintf(['\n A = ', num2str(fact)]);
    size(k2)

    k2 = k2 + k1(end,:);
    g = [g1;g2];
    k = [k1;k2];

end