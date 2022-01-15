function [fwhm] = calc_fwhm(psf)
    %-----------------------------
    % Code to calculate full width half maximum of point spread function
    %   three values for 3 orthogonal directions, x,y,z
    %   What computes here is the average radius, not actually FWHM
    %
    % Inputs:
    %   psf: NxNxN, point spread function reconstructed
    dim = size(psf);
    center_ind = floor(dim/2) + 1;
    psfz = abs(squeeze(psf(:,:,center_ind(3))));
    psfy = abs(squeeze(psf(:,center_ind(2),:)));
    psfx = abs(squeeze(psf(center_ind(1),:,:)));
    
    fwhm(1) = subfunc(psfx);
    fwhm(2) = subfunc(psfy);
    fwhm(3) = subfunc(psfz);
%     scatter3(X(:), Y(:), psfz(:)); hold on; 
%     s = scatter3(X(ptsz),Y(ptsz),psfz(ptsz),...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]);
%     s.SizeData = 100;
    
end

function [fwhm] = subfunc(psfz)
    dim = size(psfz);
    height = max(psfz(:));
    pthz = find(psfz == height,1);
    [X,Y] = meshgrid(1: dim(1), 1: dim(2));
    ptx = X(pthz); pty = Y(pthz);
    
%     figure;
%     surf(X, Y, psfz,'FaceAlpha',0.8);
%     hold on;
    
    
    
%     sf = fit([X(:), Y(:)], psfz(:),'linearinterp');
    
    [precX,precY] = meshgrid(linspace(1, dim(1), 10 * dim(1)), linspace(1, dim(2), 10 * dim(2)));
    precZ = interp2(X, Y, psfz, precX, precY);
    
    
    
%     precZ = feval(sf, [precX(:), precY(:)]);
    ptsz = find( (precZ>(0.49*height)) & (precZ<(0.51*height)) );
    fwhm = 2 * sum(sqrt((precX(ptsz) - ptx).^2 + (precY(ptsz) - pty).^2)) / length(ptsz);
%     figure;
%     surf(precX, precY, reshape(precZ, size(precX)),'FaceAlpha', 0.5);
%     hold on;
%     s = scatter3(precX(ptsz), precY(ptsz), precZ(ptsz),'r','Marker','x'); s.SizeData = 200;
%     
%     
%     ptsz = find( (psfz>(0.4*height)) & (psfz<(0.6*height)) );
%     
%     [pred_x, pred_y] = meshgrid(linspace(min(X(ptsz)), max(X(ptsz)), 100),...
%                                 linspace(min(Y(ptsz)), max(Y(ptsz)), 100));
%     
%     pred_z = feval(sf, [pred_x(:), pred_y(:)]);
%     ptsz = find( (pred_z>(0.49*height)) & (pred_z<(0.51*height)) );
%     fwhm = sum(sqrt((pred_x(ptsz) - ptx).^2 + (pred_y(ptsz) - pty).^2)) / length(ptsz);
end