function [Grs, Azi, Polar, GRCounter] = gen_Grs(Nsegs, NPhases, Nshots, ttype)
    % generate readout gradient directions after rotation
    
    GRCounter = reshape(0:(Nsegs*NPhases*Nshots-1), [], 1);
    GRCounter = reshape(GRCounter, Nshots, Nsegs, NPhases); 
    GRCounter = permute(GRCounter, [2,3,1]);    % Nsegs x Nphases x Nshots
    
    GRCounter = GRCounter(:);
    [Azi, Polar] = GoldenMeans3D(GRCounter,ttype);
            
    Grs = [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
end