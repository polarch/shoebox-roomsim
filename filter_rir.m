function rir_full = filter_rir(rir, f_center, fs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FILTER_RIR.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nBands = size(rir,2);

    if length(f_center) ~= nBands
        error('The number of bands should match the number of columns in the 2nd dimension of rir')
    end

    if nBands == 1
        rir_full = rir;
    else     
        % order of filters
        order = 1000;
        filters = zeros(order+1, nBands);
        for i=1:nBands
            if i == 1
                fl = 30;
                fh = sqrt(f_center(i)*f_center(i+1));
                wl = fl/(fs/2);
                wh = fh/(fs/2);
                w = [wl wh];
                filters(:,i) = fir1(order, w, 'bandpass');
            elseif i == nBands
                fl = sqrt(f_center(i)*f_center(i-1));
                w = fl/(fs/2);
                filters(:,i) = fir1(order, w, 'high');
            else
                fl = sqrt(f_center(i)*f_center(i-1));
                fh = sqrt(f_center(i)*f_center(i+1));
                wl = fl/(fs/2);
                wh = fh/(fs/2);
                w = [wl wh];
                filters(:,i) = fir1(order,w,'bandpass');
            end
        end

        temp_rir = [rir; zeros(order, nBands)];
        rir_filt = fftfilt(filters, temp_rir);
        rir_full = sum(rir_filt, 2);
    end
end
