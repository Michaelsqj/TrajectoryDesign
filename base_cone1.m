function [optTime, optgA, optkA, optp, optA, optArclength, undersample] = base_cone1(mat, fov, Ts, gmax, smax, readout_time, cone_angle)

    %-------------------------------------------------------------
    % Code to make a cone for Radial Cones Trajectory, Johnson MRM2016
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
    gmax        = min(gmax, (bw_readout) / (gamma * fov / 10)); %Max G/cm

    min_cost = Inf;
    %--------------------------------------------------------------
    %          find ideal p
    %--------------------------------------------------------------
    p = 0.1 : 0.3 : 5;
    t = linspace(0, 1, 1000); % t is the 'alpha' in eq[2]
    cost = zeros(numel(p), 1);

    for p_pos = 1:numel(p) % iterate to find p
        fact = 0.0;

        for step = [1 0.5 0.1 0.05 0.01 0.005]
            rtime = 0;

            while rtime < readout_time
                fact = fact + step;
                
                kmax_rad = kmax * cos(cone_angle / 180 * pi);
                k(:, 3) = t * kmax_rad;

                kTr = t.^2 * kmax * sin(cone_angle / 180 * pi);
                phi = 2 * pi * fact * t.^p(p_pos);
                k(:, 1) = cos(phi) .* kTr;
                k(:, 2) = sin(phi) .* kTr;

                [C, time, gA, s, kA, phi, sta, stb] = minTimeGradient(k, 0, 0, 0, gmax, smax, Ts);

                rtime = max(time);
            end

            fact = fact - step;
            fprintf(['\n A = ', num2str(fact), ' p = ', num2str(p(p_pos)), '\n']);
        end

        %------------Calculate Cost--------------------------------------------
        [cost(p_pos), arc_length] = costfn2(kA);
        disp(['p=', num2str(p(p_pos)), ' cost = ', num2str(cost(p_pos))]);
        
        if cost(p_pos) < min_cost
            min_cost = cost(p_pos);
            optp = p(p_pos);
            optA = fact;
            optTime = time;
            optArclength = sum(arc_length(:));
            optkA = kA;
            optgA = gA;
        end

    end
    
    undersample = optArclength / kmax;
    disp('-----------------------');
    disp(['optimal cost = ', num2str(min_cost)]);
    disp(['optimal p = ', num2str(optp)]);
    disp(['optimal A = ', num2str(optA)]);
    disp(['optimal time = ', num2str(optTime),'[ms]']);
    disp(['optimal arc length = ', num2str(optArclength), '[1/cm]']);
    disp(['undersample factor = ', num2str(undersample)]);
    disp('-----------------------');

end


function [cost, arc_length] = costfn1(kA)
    %-------------------------------------------------------------
    % Code to compute cost function as Johnson
    %
    k_r = sqrt(kA(:, 1).^2 + kA(:, 2).^2 + kA(:, 3).^2);
    rr = linspace(min(k_r), max(k_r), 1000)';

    kR(:, 1) = interp1(k_r, kA(:, 1), rr, 'linear');
    kR(:, 2) = interp1(k_r, kA(:, 2), rr, 'linear');
    kR(:, 3) = interp1(k_r, kA(:, 3), rr, 'linear');
    arc_length = sqrt(sum(diff(kR).^2, 2));

    A = rr(1:end - 1).^2; % A is sum kr.^2 in the eq[3]
    b = linsolve(A, arc_length); % b is the 'S' in the paper, scalar
    f = A * b;
    
    cost = sum((f - arc_length).^2);
end

function [cost, arc_length] = costfn2(kA)
    %-------------------------------------------------------------
    % Code to compute cost function to force the arc_length proportional
    %   to per bin weighted |kr|
    %
    k_r = sqrt(kA(:, 1).^2 + kA(:, 2).^2 + kA(:, 3).^2);
    rr = linspace(min(k_r), max(k_r), 1000)';

    kR(:, 1) = interp1(k_r, kA(:, 1), rr, 'linear');
    kR(:, 2) = interp1(k_r, kA(:, 2), rr, 'linear');
    kR(:, 3) = interp1(k_r, kA(:, 3), rr, 'linear');
    arc_length = sqrt(sum(diff(kR).^2, 2));

    A = rr(1:end - 1).^2; % A is sum kr.^2 in the eq[3]
    A = A .* gen_weights(length(A));
    b = linsolve(A, arc_length); % b is the 'S' in the paper, scalar
    f = A * b;
    
    cost = sum((f - arc_length).^2);
    
    
end

function [w] = gen_weights(len)
    %-------------------------------------------------------------
    % Code to generate weights for each point of kr^2
    % Outputs:
    %   w: value in (0,1]
    w = zeros(len,1);
    w(1:int16(len/2))   = 1;
    w(int16(len/2):end) = 1;
end
