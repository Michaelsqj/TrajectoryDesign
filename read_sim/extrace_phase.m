clear
clc
% to extract the gradient from DSV file
fpath = "";
% "ig_36x6_100hz" "radial_36x6_100hz" "rot_rad_36x6_100hz"
% "half_rad_36x6_100hz"
sim_name = "full_rad_correct_phase";
gamma = 0.4258;

% [ dsvStruct ] = dsv_readFolder( ['/Users/michael/Documents/VirtualBoxVMs/IDEAshared/sim/', sim_name] );
% XG = dsvStruct.GX; YG = dsvStruct.GY; ZG =  dsvStruct.GZ; ADC = dsvStruct.ADC; 

ADC = uncompress(fpath, sim_name, "ADC");
NC1 = uncompress(fpath, sim_name, "NC1");

ADC=ADC(1:10:end);
NC1=NC1(1:10:end);
len = min(length(NC1), length(ADC));
pts = (ADC==1); 

phase = NC1(pts);

save(sim_name+"_NC1.mat",'phase');

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