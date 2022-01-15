function [g, k, len] = minRampDown(gend, kend, smax, Ts)
    %-------------------------
    % calculate fastest ramping down gradient
    % Inputs
    %   gend: final gradient [3x1]
    %   kend: final k-space position [3x1]
    %   smax - Slew rate (G/cm/s) !!
    %   Ts -  Sampling period (s) !!
    % Outputs
    %   g:
    %   k:
    %   len:

    gzend = abs(gend(3));
    gxyend = sqrt(gend(1)^2+gend(2)^2);
    
    fprintf(['\n', 'gxyend = ', num2str(gxyend), ' gzend = ', num2str(gzend), ' smaxTs = ', num2str(smax*Ts), '\n']);
    
    gzn = ceil(gzend/(smax*Ts));
    gxyn = ceil(gxyend/(smax*Ts));
    
    len = max( gxyn, gzn );

    if (len == 0)
        g=[];
        k=[];
        return;
    elseif (len == 1)
        g = [0,0,0];
        k = [0,0,0];
        return;
    end

    temp = zeros(len+1, 3);

    temp(1: (gzn+1), 3) = linspace(gzend, 0, (gzn+1) );

    temp(:,2) = linspace(gend(:,2), 0, len+1 );
    temp(:,1) = linspace(gend(:,1), 0, len+1 );
    
    g = temp(2:end, :);

    gamma = 42.58e2;
    k = cumsum(g * Ts * gamma, 1) + kend;

end