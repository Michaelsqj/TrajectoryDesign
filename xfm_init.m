function E = xfm_init(mat, k_traj)
    dims = [mat, mat, mat, 1];
    tmpk(:,1,:) = k_traj;
    tic
    E = xfm_NUFFT(dims, [], [], tmpk);
    toc
end