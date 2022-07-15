%% load moments/ gradient from simulation
Nsegs   = 24;
Nshots  = 5;
NPhases = 3;
Navgs   = 2;

fpath = "./";
sim_name = "2022_06_16_11_11_35";

ADC = uncompress(fpath, sim_name, "ADC");
ADC = ADC(1:10:end);
XG  = uncompress(fpath, sim_name, "GRX");
YG  = uncompress(fpath, sim_name, "GRY");
ZG  = uncompress(fpath, sim_name, "GRZ");

len = min(length(ZG), length(ADC));
ADC = ADC(1:len); 
XG=XG(1:len); YG=YG(1:len); ZG=ZG(1:len);
pts = (ADC==1); 
Gx = XG(pts); Gy = YG(pts); Gz=ZG(pts);
clear XG YG ZG ADC pts
Gx  = reshape(Gx, [], Nsegs, NPhases, Navgs, Nshots) / 10;
Gy  = reshape(Gy, [], Nsegs, NPhases, Navgs, Nshots) / 10;
Gz  = reshape(Gz, [], Nsegs, NPhases, Navgs, Nshots) / 10;

%% load gradient from text file
addpath('../')
gradname = "newjohnson_grad_176_200_9920";

bwpixel         = 100;
OS              = 8;
mat             = 176;
bw_readout      = bwpixel * mat;                  % Hz
Ts              = 1e3 / bw_readout / OS;               % ms, ADC sampling time
dead_ADC_pts    = 10;
base_g          = load(fpath + gradname);
deadpts         = ceil(dead_ADC_pts * Ts / 10e-3);
T               = 10e-3;

base_g      = [base_g(:,2), base_g(:,3), base_g(:,1)];   % Phase-Read-Slice coordinate system

base_g      = [zeros(deadpts,3); base_g];

base_k(:,1) = cumtrapz(squeeze(base_g(:,1))) .* 4.258 .* T;
base_k(:,2) = cumtrapz(squeeze(base_g(:,2))) .* 4.258 .* T;
base_k(:,3) = cumtrapz(squeeze(base_g(:,3))) .* 4.258 .* T;


[R, GrRad] = calc_grmat(Nsegs, NPhases, Nshots, 1);
R   = reshape(R, Nsegs, NPhases, Nshots, 3, 3);
%% rotate the gradient/ trajectoryNsegs   = 3;
seg     = 23;
phase   = 2;
shot    = 1;

rot_g = (squeeze(R(seg, phase, shot, :, :)) * base_g')'; 
rot_k = (squeeze(R(seg, phase, shot, :, :)) * base_k')'; 

gact(:,1)   = Gx(:,seg, phase, 1, shot);
gact(:,2)   = Gy(:,seg, phase, 1, shot);
gact(:,3)   = Gz(:,seg, phase, 1, shot);

rot_g   = rot_g(1:length(gact), :);
diff_g  = rot_g - gact;
%% compare
% figure;
% scatter((1:length(rot_k))*T, rot_k(:,1)); hold on;
% scatter((1:length(rot_k))*T, rot_k(:,2)); hold on;
% scatter((1:length(rot_k))*T, rot_k(:,3)); 

figure;
scatter((1:length(rot_g))*T, rot_g(:,1)); hold on;
scatter((1:length(rot_g))*T, rot_g(:,2)); hold on;
scatter((1:length(rot_g))*T, rot_g(:,3)); 

figure;
scatter((1:length(gact))*T, gact(:,1)); hold on;
scatter((1:length(gact))*T, gact(:,2)); hold on;
scatter((1:length(gact))*T, gact(:,3)); 

figure;
scatter((1:length(diff_g))*T, diff_g(:,1)); hold on;
scatter((1:length(diff_g))*T, diff_g(:,2)); hold on;
scatter((1:length(diff_g))*T, diff_g(:,3)); 