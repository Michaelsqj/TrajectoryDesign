function include_path()
    if ~ismac
        addpath('/home/fs0/qijia/MATLAB/irt/nufft')
        addpath('/home/fs0/qijia/MATLAB/irt/utilities')
        addpath('/home/fs0/qijia/MATLAB/irt/systems')
        addpath(genpath('/home/fs0/qijia/MATLAB/irt/mex'))
        addpath('/home/fs0/qijia/MATLAB/minTimeGradient/mex-interface')
        addpath('/home/fs0/qijia/MATLAB/encoding_transforms')
        addpath('/opt/fmrib/fsl/etc/matlab')
    end
end