function [skspace] = kspace_cutoff(kspace, skmax)
    % function to reduce kmax to smaller kmax (skmax)
    % kspace: Nx3
    % skspace: Mx3
    kspace = reshape(kspace, [], 3);
    kr = sqrt(sum(kspace.^2, 2));
    pts = find(kr < skmax);
    skspace = kspace(pts,:);
end