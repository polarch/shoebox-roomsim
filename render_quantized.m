function [qIR, idx_nonzero] = render_quantized(qechogram, endtime, fs, FRACTIONAL)
%RENDER_HRTF Samples the echogram using HRTFs.
%   qechogram:      the quantized echogram structure array
%   endtime:        time in secs for desired length of the impulse response
%   fs:             sampling rate
%   FRACTIONAL:     flag for using fractional delay filters (1) for non-integer
%                   sample reflection times (more accurate), or just 
%                   quantise (0) to the closest sample (faster)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RENDER_QUANTIZED.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of grid directions of the quantization
nDirs = length(qechogram);

% render echograms
L_rir = ceil(endtime*fs);
qIR = zeros(L_rir, nDirs);
idx_nonzero = [];
for nq=1:nDirs
    tempgram = qechogram(nq);
    % omit if there are no echoes in the specific one
    if tempgram.isActive
        idx_nonzero = [idx_nonzero nq];
        % number of reflections inside the time limit
        idx_trans = find(tempgram.time<endtime, 1, 'last');
        tempgram.time = tempgram.time(1:idx_trans);
        tempgram.value = tempgram.value(1:idx_trans);

        if FRACTIONAL
            qIR(:,nq) = render_rir(qechogram(nq), endtime, fs, 1);
        else
            qIR(:,nq) = render_rir(qechogram(nq), endtime, fs, 0);
        end
    end
end

end
