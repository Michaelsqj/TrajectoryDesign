function sidelobe = calc_sidelobe(psf)
    % calc sidelobe/peak on xyz direction in central plane
    dims = size(psf);
    c = 1 + floor(dims./2);
    [pks, locs] = findpeaks(psf(:, c(2), c(3)), 1:dims(1), 'SortStr','descend');
    sidelobe(1) = pks(2)/pks(1);
    [pks, locs] = findpeaks(psf(c(1), :, c(3)), 1:dims(2), 'SortStr','descend');
    sidelobe(2) = pks(2)/pks(1);
    [pks, locs] = findpeaks(squeeze(psf(c(1), c(2), :)), 1:dims(3), 'SortStr','descend');
    sidelobe(3) = pks(2)/pks(1);
end