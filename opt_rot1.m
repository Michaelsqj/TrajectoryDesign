function opt_rot1()
    %----------------------------------
    %   This code tries to implement the 
    %       cone rotation optimization algorithm 
    %           [Johnson, MRM 2017]
    %   1. Using Pipe's algo to estimate density compensation weight for each points
    %        on K-Space
    %   while ~converge
    %       1. choose the `batch` number of shot (kj) with lowest sum of density compensation weight (maximum overlapping)
    %       for `batch` num
    %           2. gridding, get the gridded density `rhoc`.
    %           3. calculate global density with kj's contribution subtracted. `rhocs` = `rhoc` - gridding(kj)
    %           4.
    %               for rotation 0 ~ 2pi
    %                   rho(i) = inverse_gridding_to_kj (`rhoc`)
    %               end
    %           5. Choose the maximum rho(i) (minimum overlapping), `Romega` is the corresponding rotation matrix
    %           around central axis
    %           6. adding the additional contribution of density compensation weight back to `rhoc`, get new `rhocs`
    %       end
    %   end
    % Input:
    %   BaseName: 
    % Output:
    %       
end