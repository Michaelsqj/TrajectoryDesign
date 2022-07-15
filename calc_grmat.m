function [R, GrRad] = calc_grmat(Nsegs, NPhases, Nshots, ttype)
    % Output
    %    R: rotation matrix, [Nsegs*NPhases*Nshots, 3, 3]
    % calculate rotation matrix in the same way as the sequence
    
    seed = 32767.0;
    theta = gen_rand(Nsegs*NPhases*Nshots) / seed * 2.0 * pi;
    theta = reshape(theta, Nsegs, Nshots, NPhases);
    theta = permute(theta, [1,3,2]);
    theta = theta(:);
    
    [GrPRS, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, ttype);
   
    [GrPRS, GsPRS, GrRad, GsRad, R] = calc_slice(GrPRS, theta);      % R [Nsegs*NPhases*Nshots, 3, 3]

    
end