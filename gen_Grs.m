function [Grs, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, ttype)
    % generate readout gradient directions after rotation
    
    GRCounter = reshape(0:(Nsegs*NPhases*Nshots-1), [], 1);
    GRCounter = reshape(GRCounter, Nshots, Nsegs, NPhases); 
    GRCounter = permute(GRCounter, [2,3,1]);    % Nsegs x Nphases x Nshots
    
    GRCounter = GRCounter(:);
    [Azi, Polar] = GoldenMeans3D(GRCounter,ttype);

    Grs = [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
end

function [Azi, Polar] = GoldenMeans3D(N, ttype)
    % type 1: distribution on half sphere
    % type 2: distribution on whole sphere
    % Define the increments, from Chan et al
    Phis = [0.465571231876768, 0.682327803828019];

    % Calculate Polar and Azimuthal angles
    % NB. apparent error in Chan paper figure here - Beta is the angle from the
    % kz axis = Polar angle in Siemens terms

    m = N(:);
    if ttype == 1
        kz = mod(m*Phis(1),1);
    else
        kz = mod(m*Phis(1),1) * 2 - 1;
    end
    
    Polar = acos(kz); 
    Azi   = mod(m*Phis(2),1) * 2 * pi;

    if ttype == 1
        % Reverse every other line if requested
        OddIdx = logical(mod(m,2));

        % Add pi to the azimuthal angle
        Azi(OddIdx) = mod( Azi(OddIdx) + pi, 2*pi);

        % Reverse kz
        Polar(OddIdx) = acos(-kz(OddIdx)); 
        
    end

end