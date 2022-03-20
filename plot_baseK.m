function plot_baseK(base_k, base_g, kmax, gmax, T)
    % plot k-space trajectory, 1/3 region k-space trajectory, gradient waveform
    r = sqrt(sum(base_k.^2, 2));
    perfK_pts = find(r<(kmax/3) );
    figure;
    
    scatter3(base_k(:,1),base_k(:,2),base_k(:,3)); xlim([-kmax,kmax]);ylim([-kmax,kmax]);zlim([-kmax,kmax]);
    figure;
    scatter3(base_k(perfK_pts,1),base_k(perfK_pts,2),base_k(perfK_pts,3)); xlim([-kmax,kmax]);ylim([-kmax,kmax]);zlim([-kmax,kmax]);

    figure;
    scatter((1:length(base_g))*T, base_g(:,1)); hold on;
    scatter((1:length(base_g))*T, base_g(:,2)); hold on;
    scatter((1:length(base_g))*T, base_g(:,3)); 
    ylim([-gmax,gmax]);
    legend('gx','gy','gz');
    xlabel('time (us)'); ylabel('gradient (G/cm)');
end