function [kspace, grad] = gen_StackCones(mat, fov, Ts, gmax, smax, readout_time, ileaves)
    %---------------------------------
    % 
    % Output
    %   kspace: [NCols, NLines, 3]
    %   grad:   [NCols, NLines, 3]

    fov = fov/10;       % [cm]
    res = 10*fov/mat;   % [mm]
    kmax = 5/res;       % [1/cm]
    Ts = Ts/1000;       % [s]
    smax = smax * 1000; % [G/cm/s]
    LEN = readout_time / (1000*Ts);
    gamma = 42.58e2;   % Hz

    addpath('./iGurney');
    %---------------------------------
    theta         = [0+1e-5 pi/2-1e-5];

    DCF = 0;
    [g_range,k_range,nint_range,lenro] = findcone(res,fov,LEN,theta,0.1,Ts,1,smax,gmax, DCF);

    time = length(k_range) * Ts;
    
    theta_lower = min(theta);
    theta_upper = max(theta);

    %---------------------------------
    
    rmpath('./iGurney');

    Phis = [0.465571231876768, 0.682327803828019];


    kspace = zeros(length(k_range), ileaves, 3);
    grad   = zeros(length(k_range), ileaves, 3);
    gratio = 1;
    if gratio
        [~,~,~,~,~,xdir,ydir,zdir] = generate_unoptimized_radial_cones( 500e3 , 2, 12, ileaves,176,200,10e-3, 2, 30);
        for ii = 0: (ileaves-1)
            kz = mod(ii*Phis(1),1) * 2 - 1;
            polar_angle = asin(kz); 
            azim_angle  = mod(ii*Phis(2),1) * 2 * pi;
            % pos = ii+1;
            % polar_angle = asin( zdir(pos)/ ( sqrt( xdir(pos)^2 + ydir(pos)^2 + zdir(pos)^2)));
            % azim_angle = atan2( ydir(pos),xdir(pos));

            g = scale_cone(g_range, polar_angle, azim_angle, theta_lower, theta_upper);
            k = scale_cone(k_range, polar_angle, azim_angle, theta_lower, theta_upper);

            kspace(:,ii+1,:) = k;
            grad(:,ii+1,:) = g;
        end
    else
        % Calculate polar angles
        N_cones = ceil(pi * kmax * fov);
        polar_angles = linspace(-pi/2, pi/2, N_cones);

        % number of interleaves in each plane
        nint_each = (pi*fov/res*cos(polar_angles)) ./ ...
                 sqrt( 1+cos(polar_angles).^2 ./ cos(theta_lower)^2 * max(0, (pi*fov/res*cos(theta_lower)/nint_range)^2-1 ) );
        nint_each = ceil(nint_each);

        disp("Full sampling requires " + num2str(sum(nint_each)) + "interleaves");

        nint_end = cumsum(nint_each);
        nint_start = 1 + [0, cumsum(nint_each(1:end-1))];
        % calculate azimuthal angles
        azim_angles = (0:(sum(nint_each)-1)) * Phis(1) * 2 * pi;

        % select mask
        sel_mask = zeros(size(azim_angles));

        % Golden ratio redordering of polar angles
        for ii=0:(ileaves-1)
            kz = mod(ii*Phis(1),1) * 2 - 1;
            p_ind = round( (asin(kz)+pi/2) / (polar_angles(2)-polar_angles(1)) );
            p_ind = min(N_cones, max(1, p_ind));
            a_ind = nint_start(p_ind);
            
            if all(sel_mask==0,'all')
                sel_mask = zeros(size(azim_angles));
            end
            while sel_mask(a_ind)==1
                a_ind = max( 1, mod(a_ind+1, length(azim_angle)+1) );
            end

            azim_angle = azim_angles(a_ind);
            polar_angle = polar_angles(p_ind);

            g = scale_cone(g_range, polar_angle, azim_angle, theta_lower, theta_upper);
            k = scale_cone(k_range, polar_angle, azim_angle, theta_lower, theta_upper);

            kspace(:,ii+1,:) = k;
            grad(:,ii+1,:) = g;
        end
    end
end


function [k] = scale_cone(k_range, polar_angle, azim_angle, theta_lower, theta_upper )
    k(:,1) = cos(polar_angle) / cos(theta_lower) * k_range(:,1);
    k(:,2) = cos(polar_angle) / cos(theta_lower) * k_range(:,2);
    k(:,3) = sin(polar_angle) / sin(theta_upper) * k_range(:,3);

    R = [cos(azim_angle)    -sin(azim_angle)    0;
         sin(azim_angle)    cos(azim_angle)     0;
         0                  0                   1];
    k = (R * k')';

end