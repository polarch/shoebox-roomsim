function rec_echograms = rec_moduleSH(echograms, sh_orders)
%REC_MODULE Apply receiver directivity to echogram
%
%   echograms: nSrc X nRec echograms structure
%   sh_orders: vector of target SH order for each receiver position. If a 
%              scalar is passed instead of a vector, then the same order
%              applies to all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REC_MODULESH.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSrc = size(echograms,1);
nRec = size(echograms,2);
nOrders = length(sh_orders);

if nOrders == 1
    sh_orders = sh_orders*ones(nRec,1);
elseif (nRec ~= nOrders)
    error('The number of echograms should equal the number of SH encoding specs.')
end

rec_echograms = echograms;
% do nothing if all orders are zeros (omnis)
if (all(sh_orders==0) == 0)
    
    for ns=1:nSrc
        for nr=1:nRec      
            
            % Get vectors from source to receiver
            rec_vecs = echograms(ns,nr).coords;
            [azi, elev] = cart2sph(rec_vecs(:,1),rec_vecs(:,2),rec_vecs(:,3));
            polar = pi/2-elev;
            
            nSH = (sh_orders(nr)+1)^2;
            sh_gains = getSH(sh_orders(nr), [azi polar], 'real');
            rec_echograms(ns,nr).value = sh_gains .* (echograms(ns,nr).value*ones(1,nSH));
        end
    end
end

end
