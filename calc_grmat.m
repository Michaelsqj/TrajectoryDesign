function [R, endpoints] = calc_grmat(ileaves, rng)
    %-------------------------------------------------------------
    % Code to calculate rotation matrix in golden ratio ordering
    % Inputs
    %   
    % Outputs:
    %   R:
    %   endpoints:
    if nargin < 2
       rng = [-1, 1]; 
    end
    phi1    = 0.4656;
    phi2    = 0.6823;
    ind     = 1: ileaves;
    alpha   = 2*pi*phi2 .*ind;
    zdir     = mod( ind * phi1, (rng(2) - rng(1)) ) + rng(1);
    xdir    = cos(alpha).*sqrt( 1-(zdir.^2));
    ydir    = sin(alpha).*sqrt( 1-(zdir.^2));
    
    endpoints = [xdir; ydir; zdir]';
        
    rand('seed',10);
    omega_all = 2*pi*rand(ileaves,1);
%     omega_all = 2*pi*zeros(ileaves,1);
    
    R = zeros(3, 3, ileaves);
    for pos = 1:ileaves

        phi   = acos( zdir(pos)/ ( sqrt( xdir(pos)^2 + ydir(pos)^2 + zdir(pos)^2)));
        theta = pi/2+ atan2( ydir(pos),xdir(pos));

        %Unit vectors which describe the cone axis r'
        ux = xdir(pos);
        uy = ydir(pos);
        uz = zdir(pos);

        % Create a rotation matrix from the angles
        %about X
        Rphi = [ 1 0 0;
                 0 cos(phi)  -sin(phi);
                 0 sin(phi)  cos(phi)];
        Rtheta = [ cos(theta) -sin(theta) 0;
                   sin(theta)   cos(theta) 0;
                   0 0 1];

        omega = omega_all(pos);

        utu = [ ux^2  ux*uy ux*uz;
                ux*uy uy^2  uy*uz;
                ux*uz uy*uz  uz^2];
        usk = [ 0 -uz uy;
                uz  0 -ux;
                -uy ux 0];
        Raxis = utu + cos(omega)*(eye(3)-utu) + sin(omega)*usk;
        Rnet = Raxis*Rtheta*Rphi;

        R(:,:,pos) = Rnet;
    end
end