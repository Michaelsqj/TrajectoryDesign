%% Gradient acoustic resonance visualiser
% Reads in simulator gradient files and outputs waveforms & FFT
% PJL 25/01/22
close all
clear
%% Folder for simulation files
sim_name = 'ig_36x6_100hz';

%% Simulation settings
tw_len_ms   = 50;          %duration of sliding time window (ms)
update_rate = 1000;        %number of points to skip in visualiser frames

%% Scanner settings
TS = 1e-5;             % Gradient raster rate (s)
SMAX = 123.07;            % Slew rate limit (T/m/s)
GMAX = 23.28;             % Max grad amplitude (mT/m)
forbidden_frq=[590,1140];       % Forbidden frequencies (TERRA)
forbidden_bw=[100,220];         % Forbidden frequency ranges (TERRA)
% forbidden_frq=[580,1120];       % Forbidden frequencies (VERIO)
% forbidden_bw=[100,220];         % Forbidden frequency ranges (VERIO)

%% Conversions (using GE unit conventions for consistency with cones code)
corrSR = 100;    % Converting slew rates T/(m*s) -> G/cm/s
corrGA = 0.1;    % Converting gradient amps from mT/m -> G/cm
SMAX_Gcms = SMAX * corrSR;            % Slew rate limit
GMAX_Gcm = GMAX * corrGA;             % Max grad amplitude

min_frq=forbidden_frq-forbidden_bw./2;
max_frq=forbidden_frq+forbidden_bw./2;

%% Read in ADC
vf_txt = split(regexp(fileread("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim//"+sim_name+"/SimulationProtocol_ADC.dsv"),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
A_mat = readmatrix("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/"+sim_name+"/SimulationProtocol_ADC.dsv", "FileType", "text");

AG=[];
AG(1)=0;
mida=1;
gida=1;

while mida<(length(A_mat))
    if A_mat(mida+1)==A_mat(mida)
        repnum = A_mat(mida+2)+2;
        AG(gida:gida+repnum-1) = A_mat(mida);
        gida = gida+repnum;
        mida = mida+3;
    else 
        AG(gida)=A_mat(mida);
        gida = gida+1;
        mida = mida+1;
    end
end

% while midx<(size(X_mat,1)-1)
%     midx=midx+1;
%     gidx=gidx+1;
%     XG(gidx)=X_mat(midx);
%     if X_mat(midx)==X_mat(midx-1)
%         XG(gidx+1:gidx+X_mat(midx+1))=XG(gidx);
%         gidx=gidx+X_mat(midx+1);
%         midx=midx+1;
%     end
% end

ADC=cumsum(AG(1:end-1))./vert_scl;
ADC=ADC(1:10:end);




%% Read in X Gradient

vf_txt = split(regexp(fileread("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim//"+sim_name+"/SimulationProtocol_GRX.dsv"),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
X_mat = readmatrix("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/"+sim_name+"/SimulationProtocol_GRX.dsv", "FileType", "text");

XG=[];
XG(1)=0;
midx=1;
gidx=1;

while midx<(length(X_mat))
    if X_mat(midx+1)==X_mat(midx)
        repnum = X_mat(midx+2)+2;
        XG(gidx:gidx+repnum-1) = X_mat(midx);
        gidx = gidx+repnum;
        midx = midx+3;
    else 
        XG(gidx)=X_mat(midx);
        gidx = gidx+1;
        midx = midx+1;
    end
end

% while midx<(size(X_mat,1)-1)
%     midx=midx+1;
%     gidx=gidx+1;
%     XG(gidx)=X_mat(midx);
%     if X_mat(midx)==X_mat(midx-1)
%         XG(gidx+1:gidx+X_mat(midx+1))=XG(gidx);
%         gidx=gidx+X_mat(midx+1);
%         midx=midx+1;
%     end
% end

XG=cumsum(XG(1:end-1))./vert_scl;

%% Read in Y Gradient

vf_txt = split(regexp(fileread("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim//"+sim_name+"/SimulationProtocol_GRY.dsv"),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
Y_mat = readmatrix("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/"+sim_name+"/SimulationProtocol_GRY.dsv", "FileType", "text");

YG=[];
YG(1)=0;
midx=1;
gidx=1;

while midx<(length(Y_mat))
    if Y_mat(midx+1)==Y_mat(midx)
        repnum = Y_mat(midx+2)+2;
        YG(gidx:gidx+repnum-1) = Y_mat(midx);
        gidx = gidx+repnum;
        midx = midx+3;
    else 
        YG(gidx)=Y_mat(midx);
        gidx = gidx+1;
        midx = midx+1;
    end
end

% while midx<(size(Y_mat,1)-1)
%     midx=midx+1;
%     gidx=gidx+1;
%     YG(gidx)=Y_mat(midx);
%     if Y_mat(midx)==Y_mat(midx-1)
%         YG(gidx+1:gidx+Y_mat(midx+1))=YG(gidx);
%         gidx=gidx+Y_mat(midx+1);
%         midx=midx+1;
%     end
% end

YG=cumsum(YG(1:end-1))./vert_scl;

%% Read in Z Gradient

vf_txt = split(regexp(fileread("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim//"+sim_name+"/SimulationProtocol_GRZ.dsv"),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
vert_scl = str2double(vf_txt{2});
Z_mat = readmatrix("/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/"+sim_name+"/SimulationProtocol_GRZ.dsv", "FileType", "text");

ZG=[];
ZG(1)=0;
midx=1;
gidx=1;

while midx<(length(Z_mat))
    if Z_mat(midx+1)==Z_mat(midx)
        repnum = Z_mat(midx+2)+2;
        ZG(gidx:gidx+repnum-1) = Z_mat(midx);
        gidx = gidx+repnum;
        midx = midx+3;
    else 
        ZG(gidx)=Z_mat(midx);
        gidx = gidx+1;
        midx = midx+1;
    end
end


% while midx<(size(Z_mat,1)-1)
%     midx=midx+1;
%     gidx=gidx+1;
%     ZG(gidx)=Z_mat(midx);
%     if Z_mat(midx)==Z_mat(midx-1)
%         ZG(gidx+1:gidx+Z_mat(midx+1))=ZG(gidx);
%         gidx=gidx+Z_mat(midx+1);
%         midx=midx+1;
%     end
% end

ZG=cumsum(ZG(1:end-1))./vert_scl;


%%
% len = min(length(ZG), length(ADC));
% ADC = ADC(1:len); XG=XG(1:len); YG=YG(1:len); ZG=ZG(1:len);
% pts = (ADC==1); Gr = XG(pts); Gp = YG(pts); Gs=ZG(pts);
% Gr = reshape(Gr, [1000,4752]); Gp = reshape(Gp, [1000,4752]); Gs = reshape(Gs, [1000,4752]);


% %%
% [ dsvStruct ] = dsv_readFolder( ['/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/', 'ig_36x6_100hz'] );
% XG = dsvStruct.GX; YG = dsvStruct.GY; ZG =  dsvStruct.GZ; ADC = dsvStruct.ADC; ADC = ADC(1:10:end);

%% Visualise the sliding time window and FFT

% myVideo = VideoWriter("./" + sim_name, 'MPEG-4'); %open video file
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% myVideo.Quality=100;
% open(myVideo);

% gif(['/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim//', sim_name, '.gif']);

tw_len_us= tw_len_ms/1000/TS;     % correct for gradient raster rate

for tw_no=1:update_rate:length(XG)-tw_len_us

    XG_int=XG(tw_no:(tw_no+tw_len_us-1));
    YG_int=YG(tw_no:(tw_no+tw_len_us-1));
    ZG_int=ZG(tw_no:(tw_no+tw_len_us-1));
    ADC_int = ADC(tw_no:(tw_no+tw_len_us-1));
    
    h=figure(1);
    
    figure(1)
    subplot(3,1,1);
    plot((0:tw_len_us-1)*TS*1000,ADC_int,"r")
    xlim([0 tw_len_us*TS*1000])
    ylim([0, 2])
    subplot(3,1,2)

    plot((0:tw_len_us-1)*TS*1000,XG_int,"g")
    hold on
    plot((0:tw_len_us-1)*TS*1000,YG_int,"b")
    plot((0:tw_len_us-1)*TS*1000,ZG_int,"r")
    plot([0,(tw_len_us-1)*TS*1000],[GMAX_Gcm,GMAX_Gcm]/corrGA,":k")
    plot([0,(tw_len_us-1)*TS*1000],[-GMAX_Gcm,-GMAX_Gcm]/corrGA,":k")
    hold off
    xlim([0 tw_len_us*TS*1000])
    ylim([-GMAX_Gcm/corrGA*1.2 GMAX_Gcm/corrGA*1.2])
    xlabel("Time (ms)","Interpreter","latex","FontSize",12)
    ylabel("Gradient amp. (mT/m)","Interpreter","latex","FontSize",12)

    set(gca,"TickLabelInterpreter","latex","FontSize",12)
    title(strcat("X-Y, Z waveforms during",{' '},num2str(tw_len_us*TS*1e3),"ms sliding window"),"Interpreter","latex","FontSize",14)
    
    
    subplot(3,1,3)
    
    for ii = 1:length(max_frq)
        stem(min_frq(ii):max_frq(ii),repmat(3,1,max_frq(ii)-min_frq(ii)+1),":","Color",[0.5,0.5,0.5],"LineWidth",0.2)
        hold on
    end
%     stem(min_frq(1):max_frq(1),repmat(3,1,max_frq(1)-min_frq(1)+1),":","Color",[0.5,0.5,0.5],"LineWidth",0.2)
%     hold on
%     stem(min_frq(2):max_frq(2),repmat(3,1,max_frq(2)-min_frq(2)+1),":","Color",[0.5,0.5,0.5],"LineWidth",0.2)

    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(XG_int*corrGA)./tw_len_us),"-g")
    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(YG_int*corrGA)./tw_len_us),"-b")
    plot((0:(tw_len_us-1))/TS/tw_len_us,abs(fft(ZG_int*corrGA)./tw_len_us),"-r")
    xlim([0 3500])
    ylim([0 0.5])
    hold off
    xlabel("Freq. (Hz)","Interpreter","latex")
    set(gcf,"color","w")
    set(gca,"TickLabelInterpreter","latex","FontSize",12)
    title(strcat("FFT across",{' '},num2str(tw_len_us*TS*1e3),"ms window"),"Interpreter","latex","FontSize",14)
    
    drawnow
    
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
%     gif;
end

% close(myVideo);