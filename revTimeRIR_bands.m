function rt60_bands = revTimeRIR_bands(H, fs, nBands)
%REVTIMERIR_bands

% lowest octave band
band_centerfreqs(1) = 125;
for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end

nRIR = size(H,2);
lRIR = size(H,1);

    if nBands == 1
        rt60_bands(nr,1) = revTimeRIR(H, fs);
    else     
        % order of filters
        order = 1000;
        filters = zeros(order+1, nBands);
        for i=1:nBands
            if i == 1
                fl = 30;
                fh = sqrt(band_centerfreqs(i)*band_centerfreqs(i+1));
                wl = fl/(fs/2);
                wh = fh/(fs/2);
                w = [wl wh];
                filters(:,i) = fir1(order, w, 'bandpass');
            elseif i == nBands
                fl = sqrt(band_centerfreqs(i)*band_centerfreqs(i-1));
                w = fl/(fs/2);
                filters(:,i) = fir1(order, w, 'high');
            else
                fl = sqrt(band_centerfreqs(i)*band_centerfreqs(i-1));
                fh = sqrt(band_centerfreqs(i)*band_centerfreqs(i+1));
                wl = fl/(fs/2);
                wh = fh/(fs/2);
                w = [wl wh];
                filters(:,i) = fir1(order,w,'bandpass');
            end
        end
        for nr=1:nRIR
            RIR_nr = [H(:,nr); zeros(order, 1)];
            RIR_bands_nr = fftfilt(filters, RIR_nr);
            RIR_bands_nr = RIR_bands_nr(order/2+(1:lRIR),:);
            rt60_bands(nr,:) = revTimeRIR(RIR_bands_nr, fs);
        end
    end
end