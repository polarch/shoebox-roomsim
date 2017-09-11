function echogram = ims_coreMtx(room, source, receiver, type, typeValue)
%IMS_CORE Calculates an echogram of a rectangular space using ISM.
%   IMS_CORE calculates an echogram of a rectangular space using
%   the Image-Source Method, for a given source and receiver. Input
%   argument room should be a structure with the following fields:
%   room-->length, width, height, absorption. room.absoprtion is a 2x3 matrix
%   with absorption coefficients (broadband) for each of the walls on the
%   respective planes [x+ y+ z+; x- y- z-].
%
%   source and receiver are structures holding the coordinates of the
%   source/receiver as: source.coord = [Sx Sy Sz]. There are plans to
%   include directivity coefficients in the source structure.
%
%   Coordinates of source/receiver are specified from the left ground corner
%   of the room:
%                ^x
%              __|__    _
%             |  |  |   |
%             |  |  |   |
%          y<----.  |   | l
%             |     |   |
%             |     |   |
%             o_____|   -
%
%             |-----|
%                w
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMS_CORE.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALISE ENVIRONMENT

% room dimensions
l = room(1);
w = room(2);
h = room(3);

% source coordinates
if any(source<0)
    error('Source coordinates should be a positive number')
end
if (source(1)>l)||(source(2)>w)||(source(3)>h)
    error('Source coordinates out of bounds')
end
% move origin to the centre of the room
src(1) = source(1) - l/2;
src(2) = w/2 - source(2);
src(3) = source(3) - h/2;

% receiver coordinates
if any(receiver<0)
    error('Receiver coordinates should be a positive number')
end
if (receiver(1)>l)||(receiver(2)>w)||(receiver(3)>h)
    error('Receiver coordinates out of bounds')
end
% move origin to the centre of the room
rec(1) = receiver(1) - l/2;
rec(2) = w/2 - receiver(2);
rec(3) = receiver(3) - h/2;

if isequal(lower(type), 'maxorder')
    maxOrder = typeValue;
    echogram = ims_coreN(room, src, rec, maxOrder);
    
elseif isequal(lower(type), 'maxtime')
    maxDelay = typeValue;
    echogram = ims_coreT(room, src, rec, maxDelay);
    
else
    error('wrong type of ims calculation, type should be either "maxOrder" or "maxTime"')
end

% sort reflections according to propagation time
[echogram.time, idx] = sort(echogram.time,'ascend');
echogram.value = echogram.value(idx);
echogram.order = echogram.order(idx, :);
echogram.coords = echogram.coords(idx, :);

end



function reflections = ims_coreN(room, src, rec, N)

    % yeah, speed of sound...
    c = 343;

    % i,j,k indices for calculation in x,y,z respectively
    [I, J, K] = ndgrid(-N:N, -N:N, -N:N);
    % vectorize
    I = I(:);
    J = J(:);
    K = K(:);
    % compute total order and select only valid incides up to order N
    s_ord = abs(I) + abs(J) + abs(K);
    I = I(s_ord<=N);
    J = J(s_ord<=N);
    K = K(s_ord<=N);

    % image source coordinates with respect to receiver
    s_x = I*room(1) + (-1).^I*src(1) - rec(1);
    s_y = J*room(2) + (-1).^J*src(2) - rec(2);
    s_z = K*room(3) + (-1).^K*src(3) - rec(3);
    % distance
    s_d = sqrt(s_x.^2 + s_y.^2 + s_z.^2);
    % reflection propagation time
    s_t = s_d/c;
    % reflection propagation attenuation - if distance is <1m
    % set at attenuation at 1 to avoid amplification
    s_att = zeros(size(s_d));
    s_att(s_d<=1) = 1;
    s_att(s_d>1) = 1./s_d(s_d>1);

    % write to echogram structure
    reflections.value = s_att;
    reflections.time = s_t;
    reflections.order = [I J K];
    reflections.coords = [s_x s_y s_z];

end

function reflections = ims_coreT(room, src, rec, maxTime)

% yeah, speed of sound...
c = 343;

% find order N that corresponds to maximum distance
d_max = maxTime*c;
Nx = ceil(d_max/room(1));
Ny = ceil(d_max/room(2));
Nz = ceil(d_max/room(3));

    % i,j,k indices for calculation in x,y,z respectively
    [I, J, K] = ndgrid(-Nx:Nx, -Ny:Ny, -Nz:Nz);
    % vectorize
    I = I(:);
    J = J(:);
    K = K(:);
    % image source coordinates with respect to receiver
    s_x = I*room(1) + (-1).^I*src(1) - rec(1);
    s_y = J*room(2) + (-1).^J*src(2) - rec(2);
    s_z = K*room(3) + (-1).^K*src(3) - rec(3);
    % distance
    s_d = sqrt(s_x.^2 + s_y.^2 + s_z.^2);

    % bypass image sources with d>dmax
    I = I(s_d<d_max);
    J = J(s_d<d_max);
    K = K(s_d<d_max);
    s_x = s_x(s_d<d_max);
    s_y = s_y(s_d<d_max);
    s_z = s_z(s_d<d_max);
    s_d = s_d(s_d<d_max);
    
    % reflection propagation time
    s_t = s_d/c;
    % reflection propagation attenuation - if distance is <1m
    % set at attenuation at 1 to avoid amplification
    s_att = zeros(size(s_d));
    s_att(s_d<=1) = 1;
    s_att(s_d>1) = 1./s_d(s_d>1);
    
    % write to echogram structure
    reflections.value = s_att;
    reflections.time = s_t;
    reflections.order = [I J K];
    reflections.coords = [s_x s_y s_z];
    
end
