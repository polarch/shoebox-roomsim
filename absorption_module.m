function abs_echograms = absorption_module(echogram, alpha, limits)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ABSORPTION_MODULE.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nbands = size(alpha,1);

% if (Necho ~= Nbands)
%     error(['The number of echograms should equal the number of absorption ' ...
%         'bands. Echograms is an array of echogram structures, '...
%         'while absorption alpha is a matrix with absorption values for the' ...
%         'six walls per row with one row per frequency band.'])
% end

if nargin<3
    abs_echograms(1:Nbands) = echogram;
else
    for nb = 1:Nbands
        idx_limit = find(echogram.time<limits(nb), 1, 'last');

        abs_echograms(nb).time = echogram.time(1:idx_limit);
        abs_echograms(nb).value = echogram.value(1:idx_limit,:);
        abs_echograms(nb).order = echogram.order(1:idx_limit, :);
        abs_echograms(nb).coords = echogram.coords(1:idx_limit, :);
    end
end

for nb = 1:Nbands
    
    % absorption coefficients for x, y, z walls per frequency
    a_x = alpha(nb, 1:2);
    a_y = alpha(nb, 3:4);
    a_z = alpha(nb, 5:6);
    % reflection coefficients
    r_x = sqrt(1 - a_x);
    r_y = sqrt(1 - a_y);
    r_z = sqrt(1 - a_z);
    
    % split 
    i = abs_echograms(nb).order(:, 1);
    j = abs_echograms(nb).order(:, 2);
    k = abs_echograms(nb).order(:, 3);
    
    i_even = i(mod(i,2)==0);
    i_odd = i(mod(i,2)~=0);
    i_odd_pos = i_odd(i_odd>0);
    i_odd_neg = i_odd(i_odd<0);
    
    j_even = j(mod(j,2)==0);
    j_odd = j(mod(j,2)~=0);
    j_odd_pos = j_odd(j_odd>0);
    j_odd_neg = j_odd(j_odd<0);
    
    k_even = k(mod(k,2)==0);
    k_odd = k(mod(k,2)~=0);
    k_odd_pos = k_odd(k_odd>0);
    k_odd_neg = k_odd(k_odd<0);
    
    % find total absorption coefficients by calculating the
    % number of hits on every surface, based on the order per dimension
    abs_x = zeros(size(abs_echograms(nb).time));
    abs_x(mod(i,2)==0) = r_x(1).^(abs(i_even)/2) .* r_x(2).^(abs(i_even)/2);
    abs_x(mod(i,2)~=0 & i>0) = r_x(1).^ceil(i_odd_pos/2) .* r_x(2).^floor(i_odd_pos/2);
    abs_x(mod(i,2)~=0 & i<0) = r_x(1).^floor(abs(i_odd_neg)/2) .* r_x(2).^ceil(abs(i_odd_neg)/2);

    abs_y = zeros(size(abs_echograms(nb).time));
    abs_y(mod(j,2)==0) = r_y(1).^(abs(j_even)/2) .* r_y(2).^(abs(j_even)/2);
    abs_y(mod(j,2)~=0 & j>0) = r_y(1).^ceil(j_odd_pos/2) .* r_y(2).^floor(j_odd_pos/2);
    abs_y(mod(j,2)~=0 & j<0) = r_y(1).^floor(abs(j_odd_neg)/2) .* r_y(2).^ceil(abs(j_odd_neg)/2);
    
    abs_z = zeros(size(abs_echograms(nb).time));
    abs_z(mod(k,2)==0) = r_z(1).^(abs(k_even)/2) .* r_z(2).^(abs(k_even)/2);
    abs_z(mod(k,2)~=0 & k>0) = r_z(1).^ceil(k_odd_pos/2) .* r_z(2).^floor(k_odd_pos/2);
    abs_z(mod(k,2)~=0 & k<0) = r_z(1).^floor(abs(k_odd_neg)/2) .* r_z(2).^ceil(abs(k_odd_neg)/2);
    
    s_abs_tot = abs_x .* abs_y .* abs_z;
    % final amplitude of reflection
    abs_echograms(nb).value = abs_echograms(nb).value .* (s_abs_tot*ones(1,size(abs_echograms(nb).value,2)));
       
end
