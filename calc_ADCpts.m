function [gA, kA] = calc_ADCpts(g, T, Ts, LEN)
    %-----------------------------------------
    % linear interpolate the gradient raster points to generate ADC sampling points
    % Inputs:
    %   g: gradient array [Nx3]
    %   T: gradient raster time, (ms)
    %   Ts: ADC sampling points, (ms)
    %   LEN: readout length
    gamma = 42.58e2;   % Hz
    N = size(g, 1);
    grad_points = (1 : N) .* T;
    read_points = (1 : LEN) .* Ts;
    gA(:,1) = interp1(grad_points, g(:,1), read_points);
    gA(:,2) = interp1(grad_points, g(:,2), read_points);
    gA(:,3) = interp1(grad_points, g(:,3), read_points);

    kA = cumsum(gA .* gamma .* Ts);
end