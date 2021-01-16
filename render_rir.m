function IR = render_rir(echogram, endtime, fs, FRACTIONAL, RAND_IMS, rand_ims_dels)
%RENDER_RIR Samples the echogram to a specified sample rate.
%   echogram:       the echogram structure
%   endtime:        time in secs for desired length of the impulse response
%   fs:             sampling rate
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RENDER_RIR.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5, RAND_IMS = 0; end
if nargin<4, FRACTIONAL = 1; end

% number of reflections inside the time limit
idx_trans = find(echogram.time<endtime, 1, 'last');

if FRACTIONAL    
    % get lagrange interpolating filter of order 100 (filter length 101)
    order = 100;
    L_frac = order + 1;
    h_offset = 50;
    h_idx = 1 - 51:L_frac - 51;
    
    % make a filter table for quick access for quantized fractional samples
    fractions = 0:0.01:1;
    H_frac = lagrange(order, 50 + fractions);    

    % initialise array
    tempIR = zeros(ceil(endtime*fs)+2*h_offset, size(echogram.value,2));

    % render impulse response
    for i=1:idx_trans

        % select appropriate fractional delay filter
        if RAND_IMS
            refl_idx = floor((echogram.time(i)+rand_ims_dels(i))*fs) + 1;
            refl_frac = mod((echogram.time(i)+rand_ims_dels(i))*fs,1);
        else 
            refl_idx = floor(echogram.time(i)*fs) + 1;
            refl_frac = mod(echogram.time(i)*fs,1);
        end
        [~, filter_idx] = min(abs(refl_frac-fractions));
        h_frac = H_frac(:, filter_idx);
%        h_frac = lagrangeK(order, 50 + refl_frac);
        tempIR(h_offset+refl_idx+h_idx,:) = tempIR(h_offset+refl_idx+h_idx,:) + ...
            h_frac * echogram.value(i,:);
    end

    IR = tempIR(h_offset+1:end - h_offset,:);
    
else
    % initialise array
    IR = zeros(ceil(endtime*fs), size(echogram.value,2));
    % quantized indices of reflections
    if RAND_IMS
        refl_idx = round((echogram.time(1:idx_trans)+rand_ims_dels(1:idx_trans))*fs) + 1;
    else
        refl_idx = round(echogram.time(1:idx_trans)*fs) + 1;
    end
    % sum reflection amplitudes for each index
    for i=1:idx_trans
        IR(refl_idx(i),:) = IR(refl_idx(i),:) + echogram.value(i,:);
    end

end
