% function main_angi(cone_type)
    close all
    clear
    clc
    include_path();
    %% define basic parameters
    % System
    gmax    = 2.328/11;                        % G/cm
    smax    = 12.307/11;                       % G/cm/ms
    T       = 10e-3;                        % gradient raster time

    % Design Parameters
    bwpixel = 100;                      % Hz, bandwidth per pixel
    mat     = 64;              % matrix size
    fov     = 200;                      % mm
    res     = fov / mat;                % mm
    kmax    = 5 / res;                  % [1/cm]

    OS      = 4;                        % oversampling factor
    bw_readout      = bwpixel * mat;                  % Hz
    Ts              = 1e3 / bw_readout / OS;               % ms, ADC sampling time
    NCols           = mat * OS;
    readout_time    = NCols * Ts;                           % ms
    grad_time       = ceil(readout_time / T) * T;
    dead_pts        = 0;            % dead dwell time before the actual gradient
    dead_time       = dead_pts * T;

%     gmax            = min(gmax, 1/(fov/10)/(4.258*Ts));
    % Sequence parameters
    NSegs           = 36;
    NShots          = 6;
    NPhases         = 1;            % for simplicity, only simulate 1 phase here
    NLines          = NSegs * NShots;
    under_factor    = NLines / (mat*mat*pi/2);               % undersampling factor

    cone_types  = {'Johnson', 'Gurney', 'improvedGurney', 'VDGurney', 'radial', 'modified_johnson','gurney_johnson','johnson_gurney','cone9'};
    cone_type   = 6;
    full_spoke  = 1;
    
    angle_factor    = 2;                % multiple times of the original calculated angles
    cone_area       = 4 * pi / NLines;              % changeable
    cone_angle      = min( angle_factor * acosd( 1 - cone_area / 2 / pi ), 20)
%     cone_angle      = 10

    %% generate the ktraj: stack of cones / rotated cones
    org_type = 1;   % organization type, rotated cones
    if org_type == 1
        [R, GrRad] = calc_grmat(NSegs, NPhases, NShots, 2);
        if full_spoke
            [~, base_g, base_k] = gen_base_cone(mat, fov, T, gmax, smax, grad_time/2, cone_angle, cone_type);
            base_g = [-base_g; base_g];
            base_k = [-base_k; base_k];
        else
            [~, base_g, base_k] = gen_base_cone(mat, fov, T, gmax, smax, grad_time - dead_time, cone_angle, cone_type);
        end
        
        base_g = [base_g(:,2), base_g(:,3), base_g(:,1)];   % Phase-Read-Slice coordinate system
        base_k = [base_k(:,2), base_k(:,3), base_k(:,1)];   % Phase-Read-Slice coordinate system
%         base_g = [zeros(dead_pts,3); base_g
        
        [base_g, base_k] = calc_ADCpts(base_g, base_k, T, Ts, NCols);
        
        for ii = 1: (NSegs*NPhases*NShots)
            kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k')';
        end
        
        kspace = reshape(kspace, NCols, NSegs, NPhases, NShots, 3);
        kspace = kspace./kmax.*pi;

        for ii = 1:NPhases
            k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
        end
    else
        % stack of cones
        [k_traj(:,1,:), len] = gen_StackCones(mat, fov, Ts, gmax, smax, readout_time, NSegs*NShots);
        k_traj = k_traj./kmax.*pi;
    end
    
    
    %% calculate SNR, SNReff
    under_factor = ceil(mat*mat*pi/2/NLines);
    [w, kr, SNReff] = calc_dens(squeeze(k_traj(:,1,:)), mat, under_factor, pi);
    densfig = plotMeanSDSlidingWindow(kr,abs(w),max(kr)/10,20,max(kr)*0.9, SNReff); 
    %%
    E = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:));
    %% recon image
    dirname = './';
    outdir  = "result/";

%     fname = 'bSSFP_pad_thresh_small';
%     load([dirname, fname, '.mat'], 'img');
%     img = imresize3(img, [mat, mat, mat]);
%     recon_img = reshape(E.mtimes2(img), mat, mat, mat);
%     save_avw(abs(recon_img), outdir+num2str(cone_type)+"_"+fname+".nii.gz",'d',[res,res,res]);

    fname = 'radial_dn';
    load([dirname, fname, '.mat'], 'img');
    img = imresize3(img, [mat, mat, mat]);
    kimg = E * img;
%     recon_img = E.iter(kimg, @pcg, 1e-4, 50, [1,1,1,0]);
%     recon_img = reshape(abs(recon_img), mat,mat,mat);
    recon_img = reshape(E.mtimes2(img), mat, mat, mat);
    save_avw(abs(recon_img), outdir+num2str(cone_type)+"_"+fname+".nii.gz",'d',[res,res,res]);

    
    %% calculate metrics
    PSF = reshape(E' * (E.w.*ones(length(k_traj),1)), [mat, mat, mat]);
    [fwhm] = calc_fwhm(abs(PSF));
    sidelobe_level = calc_sidelobe(abs(PSF));
    save_avw(abs(PSF), outdir+num2str(cone_type)+"_PSF.nii.gz",'d',[res,res,res]);

    % -----------------------------------
    %   compute SNR (pseudo replica)
    noise_recons = add_noise(E, img);
    [SNRreplica, SNRmap] = calc_SNR(1, noise_recons, img);
    [SNRaliasing] = calc_SNR(3, noise_recons, img);
    
    % -------------------
    grid = fftshift(reshape(full(E.st(1).p' * E.w(:,1)), E.st(1).Kd));
    grid = (abs(grid)-min(abs(grid(:))))./(max(abs(grid(:))) - min(abs(grid(:))));
    save_avw(abs(grid), outdir+num2str(cone_type)+"_grid_weight.nii.gz",'d',[res/2,res/2,res/2]);
    
    %%
%     fprintf('SNReff = %.3f \n', abs(SNReff) );
    fprintf('fwhm = %.2f, %.2f, %.2f\n', fwhm(1),fwhm(2),fwhm(3));
    fprintf('sidelobe_level = %.4f, %.4f, %.4f\n', sidelobe_level(1),sidelobe_level(2),sidelobe_level(3));
    fprintf('SNRreplica = %.2f\n', abs(SNRreplica));
    fprintf('effSNRreplica = %.2f\n', abs(SNRreplica/(prod(fwhm))));
    fprintf('SNRaliasing = %.2f\n', abs(SNRaliasing));
% end