function [peak_psf, mean_psf, var_psf, eng_psf] = analyze_psf(psf, nbins)
    % calculate peak psf, mean psf, var psf, total energy of psf along radius
    psf = psf./max(psf(:));
    dim = size(psf);
    [X,Y,Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
    C = floor(dim/2)+1;
    X = X-C(1);
    Y = Y-C(2);
    Z = Z-C(3);
    R = sqrt(X.^2 + Y.^2 + Z.^2);
    maxR = 0.9*sum(C.^2).^0.5;
    psf = psf(:);
    for ii = 1:nbins
        bin_min = (ii-1)/nbins*maxR;
        bin_max = ii/nbins*maxR;
        pts = find( (R>=bin_min) .* (R<=bin_max) );
        tmp = psf(pts);
        peak_psf(ii,1) = max(abs(tmp(:)));
        mean_psf(ii,1) = mean(abs(tmp(:)));
        var_psf(ii,1) = std(tmp(:));
        eng_psf(ii,1) = sum(tmp(:).^2);
    end
end