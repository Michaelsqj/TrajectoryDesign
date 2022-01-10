function [dens_hist, SNReff] = calc_dens(dcw, rdist, max_r, num_bins)
    %---------------------------------
    % Code for calculation of hist: density compensation weigth -- radius
    if nargin < 3
        max_r = 1;
    end
    if nargin < 4
        num_bins = 64;
    end
        
    
    bin_width = max_r / num_bins;
    r_bins = floor(rdist/bin_width);
    dens_hist = zeros(num_bins, 1);
    for ii = 1:num_bins
       tmp = dcw.*(r_bins==ii);
       dens_hist(ii) = sum(tmp(:)) / sum(r_bins==ii,'all');
    end
    SNReff = sum(dcw(:))/sqrt( length(dcw) * sum( dcw(:).^2 ) );
end