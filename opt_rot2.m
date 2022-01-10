function [optR, endpoints] = opt_rot2(radial_tops, base_endpoint)
    %--------------------------------------
    %  Use the electrostatic model to uniformly distribute the endpoints by 
    %   rotating each cone around the central axis.
    save('pycode/opt.mat', 'radial_tops','base_endpoint');
    system('/home/fs0/qijia/.conda/envs/pytorch-1.1-cpu_py36/bin/python3.6 ./pycode/opt.py');
    load('/home/fs0/qijia/code/SimTraj/trajectorysimulation/pycode/optR.mat','optR','endpoints');
end