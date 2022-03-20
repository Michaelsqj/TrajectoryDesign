function test_spoiler(base_k, R, spoiler, kmax)
    % Input:
    %   base_k:     Nx3
    %   R:          Nx3x3
    %   kmax:       1/cm
    
    for ii = 1:size(R,1)
        kend(ii,:) = ( squeeze(R(ii,:,:)) * base_k(end,:)' )';
        send(ii,:) = ( squeeze(R(ii,:,:)) * (base_k(end,:)'+ [0;spoiler;0]) )';
        if ii > 1
            ckend(ii,:) = csend(ii-1,:) + kend(ii,:);
            csend(ii,:) = csend(ii-1,:) + send(ii,:);
        else
            csend(ii,:) = send(ii,:);
            ckend(ii,:) = kend(ii,:);
        end
    end
    load('cmap.mat','cmap');
    figure;
    
    xlim([-kmax kmax]);
    ylim([-kmax kmax]);
    zlim([-kmax kmax]);
    for ii = 1:36
        color_id = ceil(ii / 36 * size(cmap,1));
        color_id = min(size(cmap,1), color_id);
        plot3( [0,kend(ii,1),send(ii,1)], [0,kend(ii,2),send(ii,2)], [0,kend(ii,3),send(ii,3)], 'Color',cmap(color_id,:));
        hold on;
    end
    grid on
    colorbar;
    
    figure;
    xlim([-kmax kmax]);
    ylim([-kmax kmax]);
    zlim([-kmax kmax]);
    scatter3(0,0,0,...
        'MarkerEdgeColor','r',...
        'MarkerFaceColor',[0 .75 .75]);
    hold on;
    for ii = 1:36
        color_id = ceil(ii / 36 * size(cmap,1));
        color_id = min(size(cmap,1), color_id);
        if ii > 1
            plot3( [csend(ii-1,1),ckend(ii,1),csend(ii,1)], [csend(ii-1,2),ckend(ii,2),csend(ii,2)], [csend(ii-1,3),ckend(ii,3),csend(ii,3)], 'Color',cmap(color_id,:) );
        else
            plot3( [0,ckend(ii,1),csend(ii,1)], [0,ckend(ii,2),csend(ii,2)], [0,ckend(ii,3),csend(ii,3)], 'Color',cmap(color_id,:) );
        end
        
        hold on;
    end
    grid on
    colorbar;
end