function noise_recons = add_noise(E, img)

    reptime = 15;
    snr = 5;
    
    imgk = E * img(:);
    Pimgk = sum(imgk.^2, 'all') ./ length(imgk(:));
    Pnoi = Pimgk * 10^(-snr/10);
    sigma = sqrt(Pnoi);
    noise_recons = zeros([size(img), reptime]);

    %----------------------------
    %   randomly generete the noise several times
    for ii = 1: reptime
        tic
        noisek = sigma * randn(size(imgk)) + 1i * sigma * randn(size(imgk));
        noise_recons(:,:,:,ii) = reshape(E'*(imgk + noisek), size(img));
        toc
    end
end