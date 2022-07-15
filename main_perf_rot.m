function main_perf_rot(cone_type, ind)
    include_path();
    
    % tmpvar = [15, 20, 30, 40, 50, 60, 70]

    setParam;
    fpath       = "/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/sim_phantom_img/"
    outpath     = "/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_23-4-22/perfusion/"
    
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    phantoms    = ["ASLDRO_perfusion", "ASLDRO_perfusion_rate"]

    cone_types  = {'Johnson', 'Gurney', 'improvedGurney', 'VDGurney', 'radial', 'modified_johnson','gurney_johnson','johnson_gurney','cone9', 'seiffert_spiral'};
    % cone_type = 5;
    basename    = char(cone_types(cone_type)) + "_mat" + num2str(mat) + "_bw" + num2str(bwpixel) + "_lines" + num2str(NLines) + "_angle" + num2str(ceil(cone_angle));
    disp(basename)

    NSegs = NSegs * 2;      % use 2X segments to reconstruct
    NLines  = NSegs * NShots;
    %% ------------------------------------------------
    %    generate the ktraj: stack of cones / rotated cones
    [R, GrRad] = calc_grmat(NSegs, NPhases, NShots, 2);
    [~, base_g, base_k] = gen_base_cone(mat, fov, T, gmax, smax, grad_time - dead_time, cone_angle, cone_type);
    
    base_g = [base_g(:,2), base_g(:,3), base_g(:,1)];   % Phase-Read-Slice coordinate system
    base_k = [base_k(:,2), base_k(:,3), base_k(:,1)];   % Phase-Read-Slice coordinate system
    base_g = [zeros(dead_pts,3); base_g];
    
    [base_g, base_k] = calc_ADCpts(base_g, base_k, T, Ts, NCols);
    
    for ii = 1: (NSegs*NPhases*NShots)
        kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k')';
    end
    
    kspace = reshape(kspace, NCols, NSegs, NPhases, NShots, 3);
    kspace = kspace./kmax.*pi;

    for ii = 1:NPhases
        k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
    end

    %% ----------------------------------------------------------------
    % mat size, k_traj : 1/3 kmax; NShots = 2
    kr = sum((k_traj.^2),3).^0.5;
    perf_pts = find(kr < (pi/3));
    k_traj = k_traj(perf_pts, :, :) .* 3;      % Take 1/3 of k-space
    mat = ceil(mat / 3);
    under_factor = mat*mat*pi/2 / NLines;
    res = fov / mat;
    
    % ------------------------------------------------
    %    calculate SNR, SNReff
    [w, kr, SNReff] = calc_dens(squeeze(k_traj(:,1,:)), mat, under_factor, pi);
    densfig = plotMeanSDSlidingWindow(kr,abs(w),max(kr)/10,64,max(kr)*0.9, SNReff);
    saveas(densfig, outpath + basename + "_densfig.png");
    % ------------------------------------------------
    tic
    E = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:));
    toc
    %% recon image
    for fname = phantoms
        [img, ~,~,~,~] = read_avw(fpath+fname+".nii.gz");
        % load(fpath + fname, 'img');
        img = imresize3(img(:,:,:,1), [mat, mat, mat]);
        kdata = E * img;
        % % add weighting
        % kdata = reshape(kdata, NCols, NSegs, NShots);
        % weight = reshape(linspace(2, 3, NSegs), 1, NSegs, 1); 
        % kdata = reshape(kdata .* weight, [], 1);

        recon_img = reshape( E' * kdata, mat, mat, mat);
        save_avw(abs(recon_img), outpath + basename + "_" + fname + ".nii.gz", 'd', [res, res, res]);
    end
    
    %% calculate PSF
    PSF = reshape( E' * ( E.w .* ones(length(k_traj), 1) ), [mat, mat, mat]);
    % PSF = ifftshift(ifftn(E.PSF));
    [fwhm] = calc_fwhm(abs(PSF));
    sidelobe_level = calc_sidelobe(abs(PSF));
    save_avw(abs(PSF), outpath + basename + "_PSF.nii.gz", 'd', [res, res, res]./2);

    % % -----------------------------------
    % %   compute SNR (pseudo replica)
    noise_recons = add_noise(E, img);
    [SNRreplica, SNRmap] = calc_SNR(1, noise_recons, img);
    [SNRaliasing] = calc_SNR(3, noise_recons, img);
    
    % % -------------------
    grid = fftshift(reshape(full(E.st(1).p' * E.w(:,1)), E.st(1).Kd));
    grid = (abs(grid)-min(abs(grid(:))))./(max(abs(grid(:))) - min(abs(grid(:))));
    save_avw(abs(grid), outpath + basename + "_grid_weight.nii.gz",'d',[res/2,res/2,res/2]);
    
    % %%
    fprintf('SNReff = %.3f \n', abs(SNReff) );
    fprintf('fwhm = %.2f, %.2f, %.2f\n', fwhm(1),fwhm(2),fwhm(3));
    fprintf('sidelobe_level = %.4f, %.4f, %.4f\n', sidelobe_level(1),sidelobe_level(2),sidelobe_level(3));
    fprintf('SNRreplica = %.2f\n', abs(SNRreplica));
    fprintf('effSNRreplica = %.2f\n', abs(SNRreplica/(prod(fwhm))));
    fprintf('SNRaliasing = %.2f\n', abs(SNRaliasing));
    save(outpath + basename + "_metrics.mat", 'SNReff', 'fwhm', 'sidelobe_level', 'SNRreplica', 'SNRaliasing');
end
