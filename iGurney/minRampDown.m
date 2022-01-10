function [g, k, len] = minRampDown(gend, kend, smax, Ts)
    %-------------------------
    % calculate fastest ramping down gradient
    % Inputs
    %   gend: final gradient [3x1]
    %   kend: final k-space position [3x1]
    %   smax - Slew rate (G/cm/s)
    %   Ts -  Sampling period (s)
    % Outputs
    %   g:
    %   k:
    %   len:

    gzend = gend(3);
    gxyend = sqrt(gend(1)^2+gend(2)^2);

    len = max( ceil(gxyend/(smax*Ts)), ceil(gzend/(smax*Ts)) );

    temp = zeros(len, 3);

    temp(1: ceil(gzend/(smax*Ts)),3) = linspace(gzend, 0, ceil(gzend/(smax*Ts)) );

    temp(:,2) = linspace(gend(:,2), 0, len );
    temp(:,1) = linspace(gend(:,1), 0, len );
    
    g = temp(2:end, :);

    k = zeros(len, 3);

    gamma = 42.58e2;
    k = cumsum(g * Ts * gamma, 1) + kend;

end