function [k_traj, len] = gen_StackCones(mat, fov, Ts, gmax, smax, readout_time, ileaves)
    fov = fov/10;       % [cm]
    res = 10*fov/mat;   % [mm]
    kmax = 5/res;       % [1/cm]
    Ts = Ts/1000;       % [s]
    smax = smax * 1000; % [G/cm/s]
    LEN = readout_time / (1000*Ts);
    gamma = 42.58e2;   % Hz

    addpath('./iGurney');
    %---------------------------------
    theta         = [0+1e-5 pi/2-1e-5];

    DCF = 0.5;
    [g_range,k_range,nint_range,lenro] = findcone(res,fov,LEN,theta,0.1,Ts,1,smax,gmax, DCF);

    time = length(k_range) * Ts;
    
    theta_lower = min(theta);
    theta_upper = max(theta);

    %---------------------------------

    rmpath('./iGurney');

    Phis = [0.465571231876768, 0.682327803828019];

    len = length(k_range);
    k_traj = zeros(ileaves * length(k_range), 3);
    for pos = 0: (ileaves-1)
        kz = mod(pos*Phis(1),1) * 2 - 1;
        polar_angle = asin(kz); 
        azim_angle  = mod(pos*Phis(2),1) * 2 * pi;

        k(:,1) = cos(polar_angle) / cos(theta_lower) * k_range(:,1);
        k(:,2) = cos(polar_angle) / cos(theta_lower) * k_range(:,2);
        k(:,3) = sin(polar_angle) / sin(theta_upper) * k_range(:,3);

%         k = cumsum(g * Ts * gamma, 1);
        R = [cos(azim_angle) -sin(azim_angle) 0;
           sin(azim_angle)  cos(azim_angle) 0;
           0                0                 1];
        k = (R * k')';
        k = [k(:,2) k(:,3) k(:,1)];
        k_traj( pos*len+1 : (pos+1)*len, :) = k;
    end
end