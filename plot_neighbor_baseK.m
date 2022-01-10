function plot_neighbor_baseK(base_k, R, kmax, num_neighbor)
    %------------------------------------------------------
    %   Code to plot several neighboring trajectory of the base cone 
    %       for observation of the density
    
    
    % load('angi_Johnson_168_887_50X-result-gradopt.mat', 'R', 'base_k')
    % kmax = pi;
    % load('angi_Johnson_128_1287_20X-result-rd.mat','base_k')
    ileaves = size(R, 3);
    N = size(base_k, 1);
    endpoints = zeros(ileaves, 3);
    dist = zeros(ileaves,1);
    for ii = 1: ileaves
        endpoints(ii, :) = R(:,:,ii) * [0;0;1];
        dist(ii) = sum( (endpoints(ii,:) - [0,0,1]).^2);
    end
    
    [~, ind] = sort(dist,'ascend');
    
    figure;
    for ii = 1: num_neighbor
        edpts( (ii-1)*N+1: ii*N , :) = (R(:,:,ind(ii)) * (base_k') )';
        scatter3(edpts( (ii-1)*N+1: 4: ii*N , 1), edpts( (ii-1)*N+1: 4: ii*N , 2), edpts( (ii-1)*N+1: 4: ii*N , 3) );
        xlim([-kmax,kmax]);ylim([-kmax,kmax]);zlim([-kmax,kmax]);
        drawnow;
        hold on;
        pause(1);
    end
end