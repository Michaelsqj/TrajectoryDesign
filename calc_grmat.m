function [R, GrRad] = calc_grmat(Nsegs, NPhases, Nshots, ttype)
    % Output
    %    R: rotation matrix, [Nsegs*NPhases*Nshots, 3, 3]
    % calculate rotation matrix in the same way as the sequence

%     rand('seed',10);
%     theta = load("randnum");
%     theta = theta(1:Nsegs*NPhases*Nshots) / 32767.0 *2.0*pi;
    theta = gen_rand(Nsegs*NPhases*Nshots) / 32767.0 *2.0*pi;
%     theta = zeros(Nsegs*NPhases*Nshots,1);
    [GrPRS, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, ttype);
   
    [GrPRS, GsPRS, GrRad, GsRad, R] = calc_slice(GrPRS, theta);      % R [Nsegs*NPhases*Nshots, 3, 3]
end