function [time, base_g, base_k] = gen_radial(mat, fov, smax, gmax, len, Ts)
    %----------------------------------------------------------------
    % Code to generate radial spokes given readout length for comparison
    %
    % Inputs:
    %   mat = resolution
    %   fov    = field of view in [mm]
    %   smax   = G/cm/ms
    %   gmax =   G/cm
    %   len     = readout length
    %   Ts = sampling time [ms]
    
    res = fov/mat;  % [mm]
    kmax = 5/res;   % [1/cm]
    base_k = [zeros(len,1), zeros(len,1), (linspace(-kmax,kmax,len))'];
    base_g = [0, 0, 0; diff(base_k)];
    base_s = diff(base_g);
    if max(sqrt( sum(base_g.^2, 2) )) > gmax || max(sqrt( sum(base_s.^2, 2) )) > (smax*Ts)
       warning('Readout out length is too short and cannot be achieved')
    end
    time = len * Ts;
end