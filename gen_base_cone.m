function [time, g, k] = gen_base_cone(mat, fov, Ts, gmax, smax, readout_time, cone_angle, cone_type)
    %-------------------------------------------------------------
    % Code to make a base cone
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
    
    if nargin < 8
       cone_type = 1; 
    end
    switch cone_type
        % Johnson, MRM2016
        case 1
            [time, g, k] = base_cone1(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        % Gurney, MRM2006
        case 2
            [time, g, k] = base_cone2(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        % Improved Gurney
        case 3
            [time, g, k] = base_cone3(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        % TODO
        % Variable density Gurney
        case 4
            [time, g, k] = base_cone4(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        case 5
        % generate radial spokes with equal readout length
            [time, g, k] = gen_radial(mat, fov, smax, gmax, floor(readout_time/Ts), Ts);
        case 6
        % modified johnson's cone
            [time, g, k] = base_cone6(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        case 7
        % connect johnson + gurney
            [time, g, k] = base_cone7(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        case 8
            [time, g, k] = base_cone8(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        case 9
            [time, g, k] = base_cone9(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
        case 10
            [time, g, k] = base_cone10(mat, fov, Ts, gmax, smax, readout_time, cone_angle);
            % generate radial spokes with equal readout length
%             [time, g, k] = gen_radial(mat, fov, smax, gmax, floor(readout_time/Ts), Ts, 1);
        otherwise
            warning('unexpected cone type, cone type should be within 1~4')
    end

end