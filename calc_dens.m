% function [dens_hist, SNReff] = calc_dens(dcw, rdist, max_r, num_bins)
%     %---------------------------------
%     % Code for calculation of hist: density compensation weigth -- radius
%     if nargin < 3
%         max_r = 1;
%     end
%     if nargin < 4
%         num_bins = 64;
%     end
%         
%     
%     bin_width = max_r / num_bins;
%     r_bins = floor(rdist/bin_width);
%     dens_hist = zeros(num_bins, 1);
%     for ii = 1:num_bins
%        tmp = dcw.*(r_bins==ii);
%        dens_hist(ii) = sum(tmp(:)) / sum(r_bins==ii,'all');
%     end
%     SNReff = sum(dcw(:))/sqrt( length(dcw) * sum( dcw(:).^2 ) );
% end

function [w, kr, SNReff] = calc_dens(k_traj, mat, undersample, kmax)
    MtxSz = ceil([mat, mat, mat]./sqrt(undersample));
    ksize = [6,6,6];
    Kd    = MtxSz.*2;
    nshift = ceil(MtxSz./2);
    st = nufft_init(k_traj./kmax.*pi, MtxSz, ksize, Kd, nshift);
    P = st.p;
    clear st
    w = ones([size(k_traj,1) 1]);
    for kk=1:20
        disp(['Iteration ' num2str(kk)])
        tmp = P * (P' * w);
        w = w ./ real(tmp);
    end
    
    kr = sqrt(sum(k_traj.^2,2));
    SNReff = sum(w(:))/sqrt( length(w) * sum( w(:).^2 ) );
end