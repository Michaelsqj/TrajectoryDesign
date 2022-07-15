function gen_sim_kdata()
    % generate kspace data for simulation
    include_path();
    
    cone_type=6;
    % generate k-space trajectory
    setParam;

    fpath       = "/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/sim_phantom_img/"
    outpath     = "/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_21-6-22/"
    
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
        base_kr = sum(base_k_adc.^2, 2).^0.5;
        [Azi, Polar] = GoldenMeans3D(1:length(base_kr)*NLines, 2);
        krs = repmat(reshape(base_kr,1,[]), NLines,1);
        krs = reshape(krs, [],1);
        kspace = krs .* [sin(Azi).*sin(Polar), cos(Azi).*sin(Polar), cos(Polar)];
        kspace = reshape(kspace, [],1,3);
        k_traj = kspace./kmax.*pi;
        %----------------------------

        % for ii = 1: (NSegs*NPhases*NShots)
        %     kspace(:, ii, :) = (squeeze(R(ii,:,:)) * base_k_adc')';
        % end
        % kspace = reshape(kspace, NCols, NSegs, NPhases, NShots, 3);
        % kspace = kspace./kmax.*pi;
    
        % for ii = 1:NPhases
        %     k_traj(:,ii,:) = reshape(kspace(:,:,ii,:,:), [NCols*NLines, 1, 3]);
        % end

        % k_traj = reshape(kspace_cutoff(k_traj(:,1,:), pi*cosd(cone_angle)), [], 1, 3);
        % k_traj = k_traj./cosd(cone_angle);

        tic
        E1 = xfm_NUFFT([mat, mat, mat, 1], [], [], k_traj(:,1,:), 'table', true);
        toc
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
    
    
    % NUFFT img --> kspace, ./weighting
    % load img
    q=matfile('/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_23-6-22/img.mat');
    img=q.img;
    img=imresize3(img, [mat,mat,mat]);
    % load sens maps
    q=matfile('/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_23-6-22/sens.mat');
    sens=q.sens;
    nc=size(sens,4);
    for ii=1:nc
        s(:,:,:,ii)=imresize3(sens(:,:,:,ii), [mat,mat,mat]);
        tmp=img.*s(:,:,:,ii);
        kdata(:,ii)= E1*tmp;
        % grid
        kgrid(:,:,:,ii) = reshape(feval(E1.st(1).interp_table_adj, E1.st(1), E1.w(:,1).*kdata(:,ii)), E1.st(1).Kd);
    end
    
    % save out
    q=matfile('/home/fs0/qijia/scratch/exp_data/TrajectoryDesign/exp_23-6-22/kgrid');
    q.img=img;
    q.sens=s;
    q.kgrid=kgrid;
    q.kspace=k_traj;
end