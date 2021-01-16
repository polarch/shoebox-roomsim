function rec_echograms = rec_moduleMic(echograms, mic_specs)
%REC_MODULE Apply receiver directivity to echogram
%
%   echograms: nSrc X nRec echograms structure
%   mic_specs: [x_1 y_1 z_1 a_1; ...; x_nRec y_nRec z_nRec a_nRec];
%               x,y,z unit vector showing the orientation of each mic
%               a=0~1 is the directivity coeff d(theta) = a + (1-a)*cos(theta)
%               a=1 is omni, a=0 is dipole, a=0.5 is cardioid, etc...
%               you can also leave empty for omni receivers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REC_MODULE.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSrc = size(echograms,1);
nRec = size(echograms,2);
if nargin<2 || isempty(mic_specs)
    mic_specs = ones(nRec,4);
end
nSpecs = size(mic_specs,1);

if (nRec ~= nSpecs)
    error(['The number of echograms should equal the number of receiver directivities.'])
end

mic_vecs = mic_specs(:,1:3);
mic_coeffs = mic_specs(:,4);

rec_echograms = echograms;
% do nothing if all receivers are omnis
if ~all(mic_coeffs==1)
    
    for ns=1:nSrc
        for nr=1:nRec      
            
            nRefl = length(echograms(ns,nr).value);
            % Get vectors from source to receiver
            rec_vecs = echograms(ns,nr).coords;
            rec_vecs = rec_vecs./(sqrt(sum(rec_vecs.^2,2))*ones(1,3));
            
            mic_gains = mic_coeffs(nr) + (1-mic_coeffs(nr)).*dot( rec_vecs, repmat(mic_vecs(nr,:),[nRefl 1]), 2);         
            rec_echograms(ns,nr).value = echograms(ns,nr).value .* mic_gains;
        end
    end
end

end
