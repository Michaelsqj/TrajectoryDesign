function [time, gA, kA, p, fact] = base_cone6(mat, fov, Ts, gmax, smax, readout_time, cone_angle)

    %-------------------------------------------------------------
    % Code to make a cone for Radial Cones Trajectory, Johnson MRM2016
    %   manually set p
    %
    % Inputs:
    %   mat = resolution
    %   fov    = field of view in [mm]
    %   smax   = G/cm/ms
    %   gmax =   G/cm
    %   readout_time = time for readout in ms
    %   cone_angle = angle at edge of k-space in degrees
    %   Ts = sampling time [ms]
    %
    %Outputs:
    %   k = 3d k space trajectory
    %   time = time points in s
    %   g = required gradients
    %   p_ideal = ideal p
    %   fact_ideal = ideal A
    
    res         = fov / mat; %mm
    kmax        = 10 / (2 * res);  % [1/cm]
    gamma       = 4258; %Hz/G
    bw_readout  = 1000/Ts;   % [Hz]
%     gmax        = min(gmax, (bw_readout) / (gamma * fov / 10)); %Max G/cm
    
%     readout_time = readout_time - Ts;
    
    fact = 0.0;
    p = 1;
    for step = [1 0.5 0.1 0.05 0.01 0.005]
        rtime = 0;

        while rtime < readout_time
            fact = fact + step;
            k = genkT(fact, p, kmax, cone_angle);

            [C, time, gA, s, kA, phi, sta, stb] = minTimeGradient(k, 0, 0, -1, gmax, smax, Ts);
%             [grd, krd, lenrd] = minRampDown(gA(end,:),kA(end,:), 1e3*smax,1e-3*Ts);
%             gA = [gA; grd];
%             kA = [kA; krd];
%             time = time + lenrd*Ts;
%             fprintf(['\n A = ', num2str(fact), ' p = ', num2str(p), ' time = ', num2str(time), ' length = ', num2str(size(gA,1)), ' rlen = ', num2str(lenrd), '\n']);
            rtime = max(time);
        end
        fact = fact - step;
        fprintf(['\n\n A is ',num2str(fact),'\n']);
    end
    
    k = genkT(fact, p, kmax, cone_angle);
    
    [C, time, gA, s, kA, phi, sta, stb] = minTimeGradient(k, 0, 0, -1, gmax, smax, Ts);
%     [grd, krd, lenrd] = minRampDown(gA(end,:),kA(end,:), 1e3*smax,1e-3*Ts);
%     gA = [gA; grd];
%     kA = [kA; krd];
%     time = time + lenrd*Ts;
end

function k = genkT(fact, p, kmax, cone_angle)
    t = linspace(0, 1, 1001); % t is the 'alpha' in eq[2]
    t1 = t(1:500);
    t2 = t(501:end);
    kmax_rad = kmax * cos(cone_angle / 180 * pi);
    k(:, 3) = t.^1 * kmax_rad;

    % kTr1 = t1 * 0;
    % kTr2 = (t2 - t2(1)) / t2(end) * kmax * sin(cone_angle / 180 * pi);


    kTr2 = t2.^3 * kmax * sin(cone_angle / 180 * pi);

    kTr1 = t1 ./ t2(1) * kTr2(1);

    kTr = [kTr1, kTr2];

    phi = 2 * pi * 1.0 * fact * t.^p;
    
%             phi1 = 2 * pi * 1.0 * fact * t1.^p;
%             phi2 = 2 * pi * fact * (t2-t1(end)).^p + phi1(end);
% 
%             phi = [phi1, phi2];

    k(:, 1) = cos(phi) .* kTr;
    k(:, 2) = sin(phi) .* kTr;
end