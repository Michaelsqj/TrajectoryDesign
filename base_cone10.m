function [time, gA, kA] = base_cone10(mat, fov, Ts, gmax, smax, readout_time, cone_angle)
    % implement seiffert spiral 
    % 
    
    res         = fov / mat; %mm
    kmax        = 10 / (2 * res);  % [1/cm]
    gamma       = 4258; %Hz/G
    
    m = 0.15;
    p = 1.5;
    
    max_s = 1;
    min_s = 0;
    while 1
        step = max_s;
        k = gen_seiffert(m,p,step).*kmax;
        [C, time, gA, s, kA, phi, sta, stb] = minTimeGradient(k, 0, 0, -1, gmax, smax, Ts);
        if time < readout_time
            max_s = max_s * 2;
        else
            break;
        end
    end
    
    while 1
        step = (max_s+min_s)/2;
        k = gen_seiffert(m,p,step).*kmax;
        [C, time, gA, s, kA, phi, sta, stb] = minTimeGradient(k, 0, 0, -1, gmax, smax, Ts);
        if time < readout_time
            min_s = step;
        elseif time > readout_time
            max_s = step;
        else
            break;
        end
    end
end

function [k] = gen_seiffert(m, p, s)
    t = linspace(0,1,1001);
    k = sqrt(m);
    phi = k.*s.*t;
    [rho, z, ~] = ellipj(s.*t, m.*ones(size(t)));
    x = rho.*cos(phi);
    y = rho.*sin(phi);
    r = linspace(0,1,length(t));
    dw = r.^p;
    
    k = [x(:).*dw(:) y(:).*dw(:) z(:).*dw(:)];
end