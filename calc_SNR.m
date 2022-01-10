function [SNR, SNRmap ]= calc_SNR(type, imgs, gt)
    %-------------------------------------------------------------
    % Code to calculate SNR, several methods
    %   1. add random noise for several times with same strength, and take the average output 
    %       as signal, the std as noise. (SNR for each voxel)
    %   2. Simply take the ground truth as signal, and the difference with
    %       recon image as noise
    %   3. Dilate the ground truth and use for binary mask. Take the region outside the mask
    %       as noise, and region in undilated ground truth as signal
    % Inputs
    %   type: 
    %       1. 
    %           imgs: img X C, could be stack of images
    %           gt:   ground truth
    %       2.
    %           imgs: img X C, could be stack of images
    %           gt:   ground truth
    %   
    % Outputs:
    %   SNR: P(signal) / P(noise)

    switch type
        case 1
            [SNR, SNRmap] = sub_func1(imgs, gt);
        case 2
            SNR = sub_func2(imgs);
            SNRmap = 0;
        case 3
            [SNR, mask]= sub_func3(imgs, gt);
            SNRmap = mask;
    end

end


function [SNR, SNRmap]= sub_func1(imgs, gt)
    mask_gt = (gt~=0);
    num = size(imgs, 4);
    assert(num>1, "SNR calculation method 1 requires multiple images");
    signal = sum(imgs, 4) ./ num;
    noise = std(imgs, 1, 4);

    Psig = signal.^2;
    Pnoi = noise.^2;
    Psig_gt = mask_gt.*Psig;
    Pnoi_gt = mask_gt.*Pnoi;
    SNR = 10*log10( sum(Psig_gt(:)) / sum(Pnoi_gt(:)) );
    SNRmap = 10*log10(Psig ./ Pnoi);

    % SNR = sum(mask_gt.*SNRmap,'all') ./ sum(mask_gt(:));
end


function SNR = sub_func2(imgs, gt)
    Psig = sum(gt(:).^2)./length(gt(:));
    for ii = size(imgs, 4)
        img = imgs(:,:,:,ii);
        noise = img - gt;
        Pnoi = sum(noise(:).^2)./length(noise(:));
        SNRs(ii) = Psig/Pnoi;
    end
    SNR = sum(SNRs)/length(SNRs);
end

function [SNR, mask] = sub_func3(imgs, gt)
    num = size(imgs, 4);
    % dilate the gt first 
    mask_gt = ( gt > (max(gt(:)) / 64) );
    mask_gt = repmat(mask_gt, [1,1,1,num]);
    % randomly choose 40% voxels for estimation of noise
    tmp = imdilate(gt > (max(gt(:)) / 64), strel('sphere',2));
    mask_noi = (rand(size(tmp)) < 0.4) .* (~tmp);
    mask_noi = repmat(mask_noi, [1,1,1,num]);
    
    noise = imgs.*squeeze(mask_noi);
    Pnoi = sum(noise(:).^2) / sum(mask_noi(:));
    
    signal = imgs.*squeeze(mask_gt);
    Psig = sum(signal(:).^2) / sum(mask_gt(:));
    
    mask = mask_gt(:,:,:,1) * 1 + mask_noi(:,:,:,1) * 2;

    SNR = 10*log10(Psig / Pnoi);
end