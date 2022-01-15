function [time, g, k] = base_cone3(mat, fov, Ts, gmax, smax, readout_time, cone_angle)
    %----------------------------------------------------------------
    % Code to generate base cone based on original Gurney cone, 
    % Gurney, MRM2006
    %
    % Inputs:
    %   mat = resolution
    %   fov    = field of view in [mm]
    %   smax   = G/cm/ms
    %   gmax =   G/cm
    %   readout_time = time for readout in [ms]
    %   cone_angle = angle at edge of k-space in degrees
    %   Ts = sampling time [ms]
    %
    % Outputs:
    %   time = time points in s
    %   k = 3d k space trajectory
    %   g = required gradients
    
    fov = fov/10;       % [cm]
    res = 10*fov/mat;   % [mm]
    kmax = 5/res;       % [1/cm]
    Ts = Ts/1000;       % [s]
    smax = smax * 1000; % [G/cm/s]
    cone_angle = (90 - cone_angle)/180*pi;  % azimuthal angle in radius
    LEN = readout_time / (1000*Ts);
    gamma = 42.58e2;   % Hz
%     MAXLEN = 1500;
%     NINT = 1;

%     [g,k,len] = wcc(pi/2-cone_angle,[fov fov],1,kmax,NINT,MAXLEN,Ts,[smax*Ts smax*Ts],[gmax gmax],0,[0,1]);
%     g = g(1:len,:);
%     k = k(1:len,:);
%     time = len * Ts;

    addpath('./iGurney');
    %---------------------------------
%     theta = cone_angle;
%     [g,k,nint,len] = findcone(res,fov,LEN,theta,0.01,Ts,1,smax,gmax,0, 0, [0,1], 0);
%     time = length(g) * Ts;

    %---------------------------------
    theta         = [0+1e-5 pi/2-1e-5];
    
    DCF = 0;
    [g_range,k_range,nint_range,lenro] = findcone(res,fov,LEN,theta,0.1,Ts,1,smax,gmax, DCF);

    time = length(g_range) * Ts;
    
    theta_lower = min(theta);
    theta_upper = max(theta);

    nint = (pi*fov/res*cos(cone_angle)) ./ ...
                     sqrt( 1+cos(cone_angle).^2 ./ cos(theta_lower)^2 * max(0, (pi*fov/res*cos(theta_lower)/nint_range)^2-1 ) );


    g(:,1) = cos(cone_angle) / cos(theta_lower) * g_range(:,1);
    g(:,2) = cos(cone_angle) / cos(theta_lower) * g_range(:,2);
    g(:,3) = sin(cone_angle) / sin(theta_upper) * g_range(:,3);

    % integrate the gradient with Ts to get the trajectory
    k = cumsum(g * Ts * gamma, 1);

    %---------------------------------

    rmpath('./iGurney');
end
