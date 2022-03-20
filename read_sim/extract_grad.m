clear
clc
% to extract the gradient from DSV file
fpath = "./";
% "ig_36x6_100hz" "radial_36x6_100hz" "rot_rad_36x6_100hz"
% "half_rad_36x6_100hz"
sim_name = "2022_03_13_09_01_00";

% [ dsvStruct ] = dsv_readFolder( ['/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/', sim_name] );
% XG = dsvStruct.GX; YG = dsvStruct.GY; ZG =  dsvStruct.GZ; ADC = dsvStruct.ADC; 

ADC = uncompress(fpath, sim_name, "ADC");
ADC = ADC(1:10:end);
XG  = uncompress(fpath, sim_name, "GRX");
YG  = uncompress(fpath, sim_name, "GRY");
ZG  = uncompress(fpath, sim_name, "GRZ");

sim_name = "2022_03_13_09_06_36";
ADC2 = uncompress(fpath, sim_name, "ADC");
ADC2 = ADC2(1:10:end);
XG2  = uncompress(fpath, sim_name, "GRX");
YG2  = uncompress(fpath, sim_name, "GRY");
ZG2  = uncompress(fpath, sim_name, "GRZ");


draw_time_graph(XG-XG2, YG-YG2, ZG-ZG2, ADC);

% len = min(length(ZG), length(ADC));
% ADC = ADC(1:len); XG=XG(1:len); YG=YG(1:len); ZG=ZG(1:len);
% pts = (ADC==1); Gr = XG(pts); Gp = YG(pts); Gs=ZG(pts);
% 
% grad = [Gr(:), Gp(:), Gs(:)];
% size(grad)
% 
% fid = fopen(sim_name+"_grad", 'w');
% fprintf(fid,'%f\t%f\t%f\n',grad'); 
% fclose(fid);

function XG = uncompress(fpath, sim_name, tail)
    vf_txt = split(regexp(fileread(fpath+sim_name+"/SimulationProtocol_" + tail + ".dsv"),"[^\n]*VERTFACTOR[^\r]*","match"),"=");
    vert_scl = str2double(vf_txt{2});
    X_mat = readmatrix(fpath+sim_name+"/SimulationProtocol_" + tail + ".dsv", "FileType", "text");

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

    XG=cumsum(XG(1:end-1))./vert_scl;
end

function draw_time_graph(XG, YG, ZG, ADC)
    tw_len_ms   = 50;          %duration of sliding time window (ms)
    update_rate = 1000;        %number of points to skip in visualiser frames
    TS = 1e-5;             % Gradient raster rate (s)
    SMAX = 123.07;            % Slew rate limit (T/m/s)
    GMAX = 23.28;             % Max grad amplitude (mT/m)
    forbidden_frq=[590,1140];       % Forbidden frequencies (TERRA)
    forbidden_bw=[100,220];         % Forbidden frequency ranges (TERRA)
    corrSR = 100;    % Converting slew rates T/(m*s) -> G/cm/s
    corrGA = 0.1;    % Converting gradient amps from mT/m -> G/cm
    SMAX_Gcms = SMAX * corrSR;            % Slew rate limit
    GMAX_Gcm = GMAX * corrGA;             % Max grad amplitude

    min_frq=forbidden_frq-forbidden_bw./2;
    max_frq=forbidden_frq+forbidden_bw./2;
    
    tw_len_us= tw_len_ms/1000/TS;     % correct for gradient raster rate

    for tw_no=1:update_rate:length(XG)-tw_len_us

        if (41929>tw_no && 41929<tw_no+tw_len_us-1)
            disp("")
        end
        XG_int=XG(tw_no:(tw_no+tw_len_us-1));
        YG_int=YG(tw_no:(tw_no+tw_len_us-1));
        ZG_int=ZG(tw_no:(tw_no+tw_len_us-1));
        ADC_int = ADC(tw_no:(tw_no+tw_len_us-1));

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

    end

end