%% SETUP

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
rec = [ 3.4 2.1 1.7];
nRec = size(rec,1);

% source start
xs0 = [ 7.4 2.1 1.8];
% source velocity (m/s)
s_vel = 1;
% other point in line
xs1 = [ 3.4 3.1 1.8];
% duration of source signal (5sec)
tSig = 5;
% block length for RIR interpolation (50msec)
tBlock = 0.05;
% number of blocks
nBlocks = ceil(tSig/tBlock);
% get source positions for block RIRs
nSrc = nBlocks+1;
src = zeros(nSrc,3);
for ns=1:nSrc
    src(ns,:) = xs0 + (xs1-xs0)/sqrt(sum((xs1-xs0).^2))*(ns-1)*tBlock;
end

rec_orders = [3];

%% RUN SIMULATOR %%%
%%% Echogram
tic

% limit the RIR by reflection order or by time-limit
type = 'maxTime';
maxlim = 0.5; % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (rt60(nb)<maxlim) limits(nb) = rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
% compute echograms
[abs_echograms, rec_echograms, echograms] = compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders);

%% Rendering
% In this case all the information (e.g. SH directivities) are already
% encoded in the echograms, hence they are rendered directly to discrete RIRs
fs = 48000;
sh_rirs = render_sh_rirs(abs_echograms, band_centerfreqs, fs);

toc

%% Generate sound scenes
% Each source is convolved with the respective array IRs, and summed with
% the rest of the sources to create the microphone mixed signals

[s, fs_sig] = audioread('~/Downloads/femdan_5sec.wav');
lSig = length(s);
lBuf = ceil(lSig/(nSrc-1));

nSH = (rec_orders+1)^2;
for nsh=1:nSH
    rirs_nsh = squeeze(sh_rirs(:,nsh,:,:));
    temp(:,nsh) = fftinterpartconv(s, rirs_nsh, lBuf);
end

temp = 0.5*temp/max(max(abs(temp)));
audiowrite('femdan_rev_moving_o3.wav',temp,fs)
