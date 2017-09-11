function [abs_echograms, rec_echograms, echograms] = compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders)

nRec = size(rec,1);
nSrc = size(src,1);
% limit the RIR by reflection order or by time-limit
type = 'maxTime';
% compute echogram due to pure propagation (frequency-independent)
for ns=1:nSrc
    for nr=1:nRec
        disp('')
        disp(['Compute echogram: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        % compute echogram
        echograms(ns, nr) = ims_coreMtx(room, src(ns,:), rec(nr,:), type, max(limits));
    end
end
disp('Apply SH directivites')
rec_echograms = rec_moduleSH(echograms, rec_orders);

% apply boundary absorption
for ns = 1:nSrc
    for nr=1:nRec
        disp('')
        disp(['Apply absorption: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        % compute echogram
        abs_echograms(ns, nr, :) = absorption_module(rec_echograms(ns,nr), abs_wall, limits);
    end
end
