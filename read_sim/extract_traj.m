clear
clc
% to extract the gradient from DSV file
fpath = "";
% "ig_36x6_100hz" "radial_36x6_100hz" "rot_rad_36x6_100hz"
% "half_rad_36x6_100hz"
sim_name = "ig_2022_02_19_06_40_50";
gamma = 0.4258;

% [ dsvStruct ] = dsv_readFolder( ['/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/', sim_name] );
% XG = dsvStruct.GX; YG = dsvStruct.GY; ZG =  dsvStruct.GZ; ADC = dsvStruct.ADC; 

ADC = uncompress(fpath, sim_name, "ADC");
M0X = uncompress(fpath, sim_name, "M0X");
M0Y = uncompress(fpath, sim_name, "M0Y");
M0Z = uncompress(fpath, sim_name, "M0Z");

ADC=ADC(1:10:end);
len = min(length(M0X), length(ADC));
ADC = ADC(1:len); M0X=M0X(1:len); M0Y=M0Y(1:len); M0Z=M0Z(1:len);
pts = (ADC==1); 

kx = M0X(pts); ky = M0Y(pts); kz=M0Z(pts);

k_traj = [kx(:), ky(:), kz(:)];
k_traj = k_traj.*gamma;
size(k_traj)

fid = fopen(sim_name+"_ktraj", 'w');
fprintf(fid,'%f\t%f\t%f\n',k_traj'); 
fclose(fid);

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