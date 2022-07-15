function main_angi_rot(cone_type, tempvar)
    include_path();
    
    % tmpvar = [15, 20, 30, 40, 50, 60, 70]

    setParam;

    cone_angle = tempvar;

    fpath       = "/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/sim_phantom_img/"
    outpath     = "/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_26-6-22/"
    
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    phantoms    = ["phantom_diff", "radial_dn", "bSSFP_pad_thresh", "ge_radial", "ASLDRO_perfusion"]

    cone_types  = {'Johnson', 'Gurney', 'improvedGurney', 'VDGurney', 'radial', 'modified_johnson','gurney_johnson','johnson_gurney','cone9', 'seiffert_spiral','stack_cones'};
    basename    = char(cone_types(cone_type)) + "_mat" + num2str(mat) + "_bw" + num2str(bwpixel) + "_lines" + num2str(NLines) + "_angle" + num2str(ceil(cone_angle));
    disp(basename)

    %% ------------------------------------------------
    %    generate the ktraj: stack of cones / rotated cones
    if cone_type < 11
        [R, GrRad] = calc_grmat(NSegs, NPhases, NShots, 1);
        [~, base_g, base_k] = gen_base_cone(mat, fov, T, gmax, smax, grad_time - dead_time, cone_angle, cone_type);
        
        % [base_g,~,base_k,~,R,~,~,~] = generate_unoptimized_radial_cones( bw_readout , gmax, smax, NLines,mat,fov,T, readout_time, cone_angle);
        % R = permute(R, [3,1,2]);

        base_g = [base_g(:,2), base_g(:,3), base_g(:,1)];   % Phase-Read-Slice coordinate system
        base_k = [base_k(:,2), base_k(:,3), base_k(:,1)];   % Phase-Read-Slice coordinate system
        base_g = [zeros(dead_pts,3); base_g];
        
        grad_delay = 0;
        [base_g_adc, base_k_adc] = calc_ADCpts(base_g, base_k, T, Ts, NCols, grad_delay);
        

        %----------------------------
        % base_kr = sum(base_k_adc.^2, 2).^0.5;
        % [Azi, Polar] = GoldenMeans3D(1:length(base_kr)*NLines, 2);
        % krs = repmat(reshape(base_kr,1,[]), NLines,1);
        % krs = reshape(krs, [],1);
        % kspace = krs .* [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
        % kspace = reshape(kspace, [],1,3);
        % k_traj = kspace./kmax.*pi;
        %----------------------------

        for ii = 1: (NSegs*NPhases*NShots)
            kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k_adc')';
        end
        kspace = reshape(kspace, NCols, NSegs, NPhases, NShots, 3);
        kspace = kspace./kmax.*pi;
    
        for ii = 1:NPhases
            k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
        end

        % k_traj = reshape(kspace_cutoff(k_traj(:,1,:), pi), [], 1, 3);

        tic
        E1 = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:), 'table', true);
        toc
        
        % clear kspace ktraj
        % % add gradient delay 
        % grad_delay = 10e-3;  %ms
        % [base_g_adc, base_k_adc] = calc_ADCpts(base_g, base_k, T, Ts, NCols, grad_delay);
        
        % for ii = 1: (NSegs*NPhases*NShots)
        %     kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k_adc')';
        % end
        % kspace = reshape(kspace, NCols, NSegs, NPhases, NShots, 3);
        % kspace = kspace./kmax.*pi;
    
        % for ii = 1:NPhases
        %     k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
        % end

        % tic
        % E2 = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:));
        % toc
    else
        [traj, grad] = gen_StackCones(mat, fov, T, gmax, smax, grad_time - dead_time, NSegs*NPhases*NShots);
        % traj: nlen, ileaves, 3
        disp("rot cone needs "+num2str(NSegs*NPhases*NShots) + "ileaves, stack cones needs" + num2str(size(traj,2)) + "ileaves");
        kspace = zeros([NCols, size(traj,2:3)]);
        for ii=1:size(traj,2)
           [~,kspace(:,ii,:)] = calc_ADCpts(squeeze(grad(:,ii,:)),squeeze(traj(:,ii,:)),T,Ts,NCols);
        end
        kspace = kspace./kmax.*pi;
        for ii = 1:NPhases
            k_traj(:,ii,:) = reshape(kspace, [], 1, 3);
        end
        under_factor = 1;
        tic
        E1 = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:), 'table', true);
        toc
    end
    
    % ------------------------------------------------
    %    calculate SNR, SNReff
    % under_factor = 4;
    % [w, kr, SNReff] = calc_dens(squeeze(k_traj(:,1,:)), mat, under_factor, pi);
    % densfig = plotMeanSDSlidingWindow(kr,abs(w),max(kr)/10,64,max(kr)*0.9, SNReff);
    % saveas(densfig, outpath + basename + "_densfig.png");
    % ------------------------------------------------
    % tic
    % E = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:));
    % toc
    %% recon image
    phantoms = ["MPRAGE-T1_new", "shepp-logan", "ge_radial"];
    for fname = phantoms
        if fname == "shepp-logan"
            img = phantom3d(mat);
        else
            [img, ~,~,~,~] = read_avw(fpath+fname+".nii.gz");
            img = imresize3(img(:,:,:,1), [mat, mat, mat]);
        end

        % kdata = E1 * img;
        % add weighting
        % kdata = reshape(kdata, NCols, NSegs, NShots);
        % weight = reshape(linspace(10000, 0, NSegs), 1, NSegs, 1);
        % weight = reshape(exp(-5 * (0:(NCols-1)) / NCols), NCols, 1, 1); 
        % kdata = reshape(kdata, [], 1);

        dd = reshape(E1.mtimes2(img), [], 1);
        recon_img = reshape(dd, mat, mat, mat);
        % W = wavelet3([mat, mat, mat]);
        % [recon_img] = fista(E1, dd, W, 0.1, [mat, mat, mat], 10, 1);
        % recon_img = reshape(E1.iter(dd, @pcg, 1e-4, 100, [1, 1, 1, 1]*0), mat, mat, mat);
        % recon_img = reshape(pogm_LLR(E1, dd, 0.05, [1 1 1]*5, [E1.Nd E1.Nt], 10), mat, mat, mat);
%         lambda = 1e4;
%         recon_img = reshape( fgp_xTV(E1, dd, lambda), mat, mat, mat);

        % recon_img = reshape( E1' * kdata, mat, mat, mat);
        save_avw(abs(recon_img), outpath + basename + "_" + fname + ".nii.gz", 'd', [res, res, res]);
    end
    
    %% calculate PSF
    PSF = reshape( E1' * ( E1.w .* ones(length(k_traj), 1) ), [mat, mat, mat]);
    % PSF = ifftshift(ifftn(E1.PSF));
    % [fwhm] = calc_fwhm(abs(PSF));
    % sidelobe_level = calc_sidelobe(abs(PSF));
    save_avw(abs(PSF), outpath + basename + "_PSFabs.nii.gz", 'd', [res, res, res]);
    save_avw(real(PSF), outpath + basename + "_PSFreal.nii.gz", 'd', [res, res, res]);
    save_avw(imag(PSF), outpath + basename + "_PSFimag.nii.gz", 'd', [res, res, res]);
    save_avw(angle(PSF), outpath + basename + "_PSFphase.nii.gz", 'd', [res, res, res]);

    % -----------------------------------
    %   compute SNR (pseudo replica)
%     noise_recons = add_noise(E, img);
%     [SNRreplica, SNRmap] = calc_SNR(1, noise_recons, img);
%     [SNRaliasing] = calc_SNR(3, noise_recons, img);
    
    % -------------------
    % for ii=1:NLines
    %     ind=(NCols*(ii-1)+1):(NCols*ii);
    %     spokes=zeros(NCols*NLines,1); spokes(ind,1)=1;
    %     grid = fftshift(reshape(full(E1.st(1).p' *spokes), E1.st(1).Kd));
    %     tmp(:,:,ii) = squeeze(abs(grid(:,129,:)));
    % end
    % save_avw(abs(tmp), outpath + basename + "_tmp.nii.gz",'d',[res/2,res/2,1]);
%     grid = fftshift(reshape(full(E1.st(1).p' * E1.w(:,1)), E1.st(1).Kd));
    grid = fftshift(reshape(feval(E1.st(1).interp_table_adj, E1.st(1), ones(size(E1.w(:,1)))), E1.st(1).Kd));
    grid = (abs(grid)-min(abs(grid(:))))./(max(abs(grid(:))) - min(abs(grid(:))));
    save_avw(abs(grid), outpath + basename + "_grid_weight.nii.gz",'d',[res/2,res/2,res/2]);

    mask = zeros(size(E1.w(:,1)));
    mask(1:length(base_k_adc),1) = 1;
    grid = fftshift(reshape(feval(E1.st(1).interp_table_adj, E1.st(1), E1.w(:,1).*mask), E1.st(1).Kd));
    grid = (abs(grid)-min(abs(grid(:))))./(max(abs(grid(:))) - min(abs(grid(:))));
    save_avw(abs(grid), outpath + basename + "_single_cone.nii.gz",'d',[res/2,res/2,res/2]);
    


    %%
    % fprintf('SNReff = %.3f \n', abs(SNReff) );
    % fprintf('fwhm = %.2f, %.2f, %.2f\n', fwhm(1),fwhm(2),fwhm(3));
    % fprintf('sidelobe_level = %.4f, %.4f, %.4f\n', sidelobe_level(1),sidelobe_level(2),sidelobe_level(3));
%     fprintf('SNRreplica = %.2f\n', abs(SNRreplica));
%     fprintf('effSNRreplica = %.2f\n', abs(SNRreplica/(prod(fwhm))));
%     fprintf('SNRaliasing = %.2f\n', abs(SNRaliasing));
%     save(outpath + basename + "_metrics.mat", 'SNReff', 'fwhm', 'sidelobe_level', 'SNRreplica', 'SNRaliasing');
end
