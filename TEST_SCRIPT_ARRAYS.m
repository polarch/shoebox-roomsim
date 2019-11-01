%% SETUP
datapath = '/home/politis/Devel/spatial_acoustic_analysis/data/';
addpath(datapath);

% room definition
room = [10.2 7.1 3.2];

% desired RT per octave band, and time to truncate the responses
rt60 = [1.0 0.8 0.7 0.6 0.5 0.4];
nBands = length(rt60);
% lowest octave band
band_centerfreqs(1) = 125;
for (nb=2:nBands) band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60);
% critical distance for the room
[~,  d_critical] = room_stats(room, abs_wall);

% receiver position
rec = [ 4.5 3.4 1.5;
        2.0 3.1 1.4;
        3.3 2.3 1.6];
nRec = size(rec,1);

% source positions
src = [ 6.2 2.0 1.8;
        5.8 5.0 1.9];
nSrc = size(src,1);

% % convert source directions from listener-centric to room-centric
% [src_coords(:,1), src_coords(:,2), src_coords(:,3)] = sph2cart(src_dirs(:,1)*pi/180, ...
%     src_dirs(:,2)*pi/180, src_r);
% src_coords(:,2) = -src_coords(:,2);
% src = ones(nSrc,1)*rec(1,:) + src_coords;
% % check sources
% for n=1:nSrc
%     if (src(n,1)>room(1))||(src(n,2)>room(2))||(src(n,3)>room(3))
%         error('Source coordinates out of room boundaries')
%     end
% end

%% DEFINE DIRECTIONAL IRs FOR RECEIVERS
fs = 48000;

% Receiver 1 example: Eigenmike SMA
disp('Simulating EM32 responses for grid')
mic_dirs_deg = [0 32 0 328 0 45 69 45 0 315 291 315 91 90 90 89 180 212 180 148 180 225 249 225 180 135 111 135 269 270 270 271;
                21 0 -21 0 58 35 0 -35 -58 -35 0 35 69 32 -31 -69 21 0 -21 0 58 35 0 -35 -58 -35 0 35 69 32 -32 -69]';
mic_dirs_rad = mic_dirs_deg*pi/180;
R = 0.042;
arrayType = 'rigid';
c = 343;
f_max = 16000;
kR_max = 2*pi*f_max*R/c;
array_order = ceil(2*kR_max);
L_resp = 1024;
% define grid to simulate responses (or that's coming from measurements directly)
grid = loadSphGrid('N040_M840_Octa.dat'); % 840 points uniformly distributed
grids{1} = grid.aziElev;
% compute array IRs
h_eigen = simulateSphArray(L_resp, mic_dirs_rad, grid.aziElev, arrayType, R, array_order, fs);
array_irs{1} = h_eigen;
clear h_eigen

% Receiver 2 example: Uniform circular array of 8 omnis, radius 10cm
disp('Simulating 8ch UCA responses for grid')
mic_dirs_deg = (0:360/8:360-360/8)';
mic_dirs_deg(:,2) = zeros(size(mic_dirs_deg));
R = 0.1;
[mic_xyz(:,1), mic_xyz(:,2), mic_xyz(:,3)] = sph2cart(mic_dirs_deg(:,1), mic_dirs_deg(:,2), R);

% impulse response parameters
L_resp = 256;
% define grid to simulate responses (same as for receiver 1)
grids{2} = grid.aziElev;
% simulate array using getArrayResponse()
fDirectivity = @(angle) 1; % response of omnidirectional microphone
h_uca = getArrayResponse(grid.xyz, mic_xyz, [], fDirectivity, L_resp, fs); % microphone orientation irrelevant in this case
array_irs{2} = h_uca;
clear h_uca

% Receiver 3 example: Measured HRTFs (MAKE SURE THAT THEY ARE THE SAME
% SAMPLERATE AS THE REST OR RESAMPLE)
disp('Loading measured HRTF responses')
load('APolitis_ownsurround2016_L512_N836_48k.mat','hrtf_dirs_deg_aziElev','hrirs')
grids{3} = hrtf_dirs_deg_aziElev*pi/180;
array_irs{3} = hrirs;
clear hrirs hrtf_dirs_deg_aziElev

%% RUN SIMULATOR %%%
%%% Echogram
tic

% limit the RIR by reflection order or by time-limit
type = 'maxTime';
maxlim = 1.5; % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (rt60(nb)<maxlim) limits(nb) = rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
% compute echograms
[abs_echograms, echograms] = compute_echograms_arrays(room, src, rec, abs_wall, limits);

%% Rendering
% In the array case, for each receiver position, a grid of directions
% should be provided, along with the corresponding multichannel IRs for the
% elements of the array at the grid directions. 
array_rirs = render_array_rirs(abs_echograms, band_centerfreqs, fs, grids, array_irs);

toc

%% Generate sound scenes
% Each source is convolved with the respective array IRs, and summed with
% the rest of the sources to create the microphone mixed signals

sourcepath = '4src_samples_voice_handclaps_fountain_piano.wav';
[src_sigs, fs_src] = audioread(sourcepath);
if (fs_src~=fs), resample(src_sigs, fs, fs_src); end % resample if it doesn't match the project's fs

array_sigs = apply_source_signals(array_rirs, src_sigs);
