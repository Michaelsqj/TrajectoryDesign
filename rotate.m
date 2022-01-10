function [k_traj] = rotate(base_k, R)
    %-------------------------------------------------------------
    % Code to rotate the base trajectory
    % Inputs
    %   R: 3x3xN
    %   base_k: Nx3
    len = size(base_k, 1);
    num = size(R, 3);
    k_traj = zeros( len * num, 3 ); 
    for ii = 1: num
        k_traj( (ii-1)*len+1: ii*len, : ) = base_k * (R(:,:,ii)'); 
    end
end