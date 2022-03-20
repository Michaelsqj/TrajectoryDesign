function plotGR(k_traj, ileaves, kmax, n)
    if nargin<4
        n=ileaves;
    end
    len = length(k_traj)/ileaves;
    load('cmap.mat','cmap');
    
    % Initialize video
    myVideo = VideoWriter('myVideoFile','MPEG-4'); %open video file
    myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
    myVideo.Quality=100;
    open(myVideo);
    
    figure;
    grid off;
  
    for ii = 1:n
        color_id = ceil(ii / n * size(cmap,1));
        color_id = min(size(cmap,1), color_id);
        plot3(k_traj((ii-1)*len+1:ii*len,1),k_traj((ii-1)*len+1:ii*len,2),k_traj((ii-1)*len+1:ii*len,3),'Color',cmap(color_id,:),'linewidth',3);
        xlim([-kmax,kmax]); ylim([-kmax,kmax]); zlim([-kmax,kmax]);
        colorbar;
        pause(0.25) %Pause and grab frame
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
        hold on
    end
    close(myVideo);
end