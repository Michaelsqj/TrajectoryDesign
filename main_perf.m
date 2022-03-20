clear
close all
clc

% include_path();

dirname = '';
outdir  = '';
if ~isfolder(outdir)
    mkdir(outdir)
end

fname = 'ASLDRO_perfusion';
load([dirname, fname, '.mat'], 'img', 'scales');
%----------------------------
%    Define the parameters

% System
gmax    = 2.328;                        % G/cm
smax    = 12.307;                       % G/cm/ms
T       = 10e-3;                        % gradient raster time, ms

% Design Parameters
bwpixel = 100;                      %Hz, bandwidth per pixel
mat     = size(img,1) * 3;              % matrix size
fov     = 200;                      % mm
res     = fov / mat;                % mm
kmax    = 5 / res;                  % [1/cm]

OS      = 4;                        % oversampling factor
bw_readout      = bwpixel * mat;                  % Hz
Ts              = 1e3 / bw_readout / OS;               % ms, ADC sampling time
LEN             = mat * OS;
readout_time    = LEN * Ts;                           % ms

grad_time       = ceil(readout_time / T) * T;

cone_types  = {'Johnson', 'Gurney', 'improvedGurney', 'VDGurney', 'radial', 'modified_johnson','gurney_johnson','johnson_gurney'};
% optims      = {'NOopt', 'Endopt', 'Jopt'};

angle_factor    = 5;                % multiple times of the original calculated angles
under_factor    = 50;               % undersampling factor

ileaves         = 2 * ceil(mat^2*pi/2/under_factor);                   % Number of interleaves
num_gold_recon  = ileaves;            % number of interleaves used in golden angle recon
cone_area       = 4 * pi / (num_gold_recon/2);              % changeable
cone_angle      = min( angle_factor * acosd( 1 - cone_area / 2 / pi ), 20)

fid = fopen('./logs/main_perf.txt','w');
fprintf(fid, '\n\n\n\n\n%s \n\n',datestr(now, 'dd/mm/yy-HH:MM') );
fprintf(fid, '\n%s\n\n',fname);

for cone_type = [6]

    basename = [outdir, '/',char(cone_types(cone_type)),'_',num2str(mat), '_', num2str(num_gold_recon), '_', num2str(angle_factor),'Xangle_', num2str(under_factor),'X', num2str(OS),'OS-']

    if cone_type == 5
        [R, endpoints] = calc_grmat(ileaves, [0, 1]);
    else
        [R, endpoints] = calc_grmat(ileaves, [-1, 1]);               % golden ratio ordering rotation matrix
    end

    %-----------------------------------
    %   generate the base cones
    [time, base_g, base_k] = gen_base_cone(mat, fov, T, gmax, smax, grad_time, cone_angle, cone_type);

    % interpolate the gradient raster points
    [base_g, base_k] = calc_ADCpts(base_g, base_k, T, Ts, LEN);

    %-----------------------------------
    %   rotate the base cone for 
    %     full k-space trajectory
    k_traj = rotate( base_k, R(:,:, 1: num_gold_recon) );

    % -----------------------------------
    %   perfusion, 1/3 in the center
    k_r = sqrt(sum((k_traj.^2),2));
    perf_pts = find(k_r < (kmax/3));
    k_traj_perf = k_traj(perf_pts, :);      % Take 1/3 of k-space

    %----------------------------------------------------
    %   recon PSF
%     Epsf = xfm_init(mat, k_traj_perf ./ kmax .*pi);
%     PSF = reshape(Epsf' * ones(size(k_traj_perf, 1), 1), [mat, mat, mat]);
%     [fwhm] = calc_fwhm(abs(PSF));
%     sidelobe_level = calc_sidelobe(abs(PSF));
    
    clear Epsf
    %-----------------------------------
    %   recon
    E = xfm_init(mat/3, k_traj_perf ./ (kmax/3) .* pi);
    recon_img = reshape(E' * (E * img(:)), [mat/3, mat/3, mat/3]);
    [dens_hist, SNReff] = calc_dens(E.w, sqrt(sum(k_traj_perf.^2, 2)), 1/3*kmax, 64);

    % -----------------------------------
    %   compute SNR (pseudo replica)
    noise_recons = add_noise(E, img);
    [SNRreplica, SNRmap] = calc_SNR(1, noise_recons, img);
    [SNRaliasing] = calc_SNR(3, noise_recons, img);

    clear E noise_recons

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid, '%s\n\n', basename);
    fprintf(fid, 'perf.SNReff = %.3f \n', abs(SNReff) );
    fprintf(fid, 'perf.fwhm = %.2f, %.2f, %.2f\n', fwhm(1),fwhm(2),fwhm(3));
    fprintf(fid, 'perf.sidelobe_level = %.4f, %.4f, %.4f\n', sidelobe_level(1),sidelobe_level(2),sidelobe_level(3));
    fprintf(fid, 'perf.SNRreplica = %.2f\n', abs(SNRreplica));
    fprintf(fid, 'perf.effSNRreplica = %.2f\n', abs(SNRreplica/(prod(fwhm))));
    fprintf(fid, 'perf.SNRaliasing = %.2f\n', abs(SNRaliasing));

    save_avw(abs(SNRmap), [basename, 'perf_SNRmap.nii.gz'], 'd', scales);
    save_avw(abs(PSF), [basename, 'perf_PSF.nii.gz'],'d',scales);
    save_avw(abs(recon_img), [basename, 'perf_recon.nii.gz'],'d',scales);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

fclose(fid);