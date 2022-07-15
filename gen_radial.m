function [time, base_g, base_k] = gen_radial(mat, fov, smax, gmax, len, Ts, ttype)
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
    gamma = 42.58e2;   % Hz/G
    res = fov/mat;  % [mm]
    kmax = 5/res;   % [1/cm]
    if nargin < 7
        base_k = [zeros(len,1), zeros(len,1), (linspace(-kmax,kmax,len))'];
    else
        base_k = [zeros(len,1), zeros(len,1), (linspace(0,kmax,len))'];
    end
    base_g = ones(size(base_k)) .* (base_k(2,:) - base_k(1,:)) ./ Ts ./ gamma .* 1e3;
%     base_s = diff(base_g);
%     if max(sqrt( sum(base_g.^2, 2) )) > gmax || max(sqrt( sum(base_s.^2, 2) )) > (smax*Ts)
%        warning('Readout out length is too short and cannot be achieved')
%     end
    time = len * Ts;
end