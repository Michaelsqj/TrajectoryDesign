%% define basic parameters
% System
gamma   = 4.258;            % [1/G/ms]
gfact   = 1;               % factor to limit the gradient and slew rate             
gmax    = 2.328/gfact;                        % G/cm
smax    = 12.307/gfact;                       % G/cm/ms
T       = 10e-3;                        % gradient raster time

% Design Parameters
bwpixel = 200;                      % Hz, bandwidth per pixel
mat     = 176;              % matrix size
fov     = 200;                      % mm
res     = fov / mat;                % mm
kmax    = 5 / res;                  % [1/cm]

OS      = 8;                        % oversampling factor
bw_readout      = bwpixel * mat;                  % Hz
Ts              = 1e3 / bw_readout / OS;               % ms, ADC sampling time
NCols           = mat * OS;
readout_time    = NCols * Ts;                           % ms
grad_time       = ceil(readout_time / T) * T;
dead_ADC_pts    = 0;            % dead ADC points before the actual gradient
dead_time       = T * ceil(dead_ADC_pts * Ts / T);
dead_pts        = ceil(dead_ADC_pts * Ts / T)

gmax            = min(gmax, 1/(fov/10)/(4.258*Ts))
% Sequence parameters
NSegs           = 36;
NShots          = 48;
NPhases         = 1;            % for simplicity, only simulate 1 phase here
NLines          = NSegs * NShots;
under_factor    = mat*mat*pi/2 / NLines;               % undersampling factor

angle_factor    = 5;                % multiple times of the original calculated angles
cone_area       = 4 * pi / NLines;              % changeable
cone_angle      = min( angle_factor * acosd( 1 - cone_area / 2 / pi ), 20)