% function [gA, kA] = calc_ADCpts(g, T, Ts, LEN)
%     %-----------------------------------------
%     % linear interpolate the gradient raster points to generate ADC sampling points
%     % Inputs:
%     %   g: gradient array [Nx3]
%     %   T: gradient raster time, (ms)
%     %   Ts: ADC sampling points, (ms)
%     %   LEN: readout length
%     gamma = 42.58e2;   % 1/G/s
%     N = size(g, 1);
%     T = T*1e3;  % us
%     Ts = Ts*1e3;    %us
%     grad_points = (0: (N-1)) .* T;
%     read_points = (0: (LEN-1)) .* Ts;
%     assert(grad_points(end) >= read_points(end));
%     T2 = 2;     %us
%     grad_points2 = 0: T2: grad_points(end);
%     % upsample the gradient first, 2us space
%     if T>T2
%         gA(:,1) = interp1(grad_points, g(:,1), grad_points2, 'linear');
%         gA(:,2) = interp1(grad_points, g(:,2), grad_points2, 'linear');
%         gA(:,3) = interp1(grad_points, g(:,3), grad_points2, 'linear');
%     end
%     
%     % figure; scatter((0:size(g,1)-1).*T,g'); legend('gx','gy','gz'); hold on; scatter((0:size(gA,1)-1).*T2,gA','x');
%     k = cumtrapz(gA) * gamma * T2 * 1e-6;
%     kA(:,1) = interp1(grad_points2, k(:,1), read_points, 'linear', 'extrap');
%     kA(:,2) = interp1(grad_points2, k(:,2), read_points, 'linear', 'extrap');
%     kA(:,3) = interp1(grad_points2, k(:,3), read_points, 'linear', 'extrap');
%     
%     % figure; scatter((0:size(k,1)-1).*T2,k'); legend('kx','ky','kz'); hold on; scatter((0:size(kA,1)-1).*Ts,kA','x');
% end

function [gA, kA] = calc_ADCpts(g, k, T, Ts, LEN)
    %-----------------------------------------
    % linear interpolate the gradient raster points to generate ADC sampling points
    % Inputs:
    %   g: gradient array [Nx3]
    %   T: gradient raster time, (ms)
    %   Ts: ADC sampling points, (ms)
    %   LEN: readout length
    gamma = 42.58e2;   % Hz/G
    N = size(g, 1);
    grad_points = (0: (N-1)) .* T;
    read_points = (0: (LEN-1)) .* Ts;
    gA(:,1) = interp1(grad_points, g(:,1), read_points, 'linear', 'extrap');
    gA(:,2) = interp1(grad_points, g(:,2), read_points, 'linear', 'extrap');
    gA(:,3) = interp1(grad_points, g(:,3), read_points, 'linear', 'extrap');
    
    kA(:,1) = interp1(grad_points, k(:,1), read_points, 'linear', 'extrap');
    kA(:,2) = interp1(grad_points, k(:,2), read_points, 'linear', 'extrap');
    kA(:,3) = interp1(grad_points, k(:,3), read_points, 'linear', 'extrap');


%     kA = cumsum(gA) .* gamma .* Ts .* 1e-3 + k(1,:);
end