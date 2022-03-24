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