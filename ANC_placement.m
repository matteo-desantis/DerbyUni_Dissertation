% ***** ANC PLACEMENT *****

% ***** Load Network *****
% Networks are stored in the network folder
% define symmetry = 1 if the target sound field is symmetric, symmetry = 0 otherwise
symmetry = 1;
if symmetry == 1
net = coder.loadDeepLearningNetwork( ...
    'Network/CNN_4D_norm_sym_train5000_5_rmse043_drop00_19noalias.mat');
elseif symmetry == 0
net = coder.loadDeepLearningNetwork( ...
    'Network/CNN_4D_norm_non_sym_train5000_5_rmse048.mat');
end

% ***** Generate source and control region grids *****
% target field sources
x_main_sources = 0;
y_main_sources = -7.0:2.0:7.0;
% control
x_control = 20.0:2:30.0;
y_control = -10.0:2:10.0;

% control for surf usage
x_control_surf = 19.0:2:31.0;
y_control_surf = -11.0:2:11.0;

% number of secondary sources
x_sources = 0.0:1.0:5.0;
y_sources = -7.5:1.0:7.5;
secondary_number = numel(x_sources)*numel(y_sources);

% ***** Main Array Grid *****
x_main = -50;
y_main = y_main_sources;
main_array = ones(numel(y_main), 1);
x_scatter3 = x_main*ones(numel(y_main), 1);
% main array gain
gain = 1;
main_array = gain*main_array;

% ***** 1/12 octave spaced frequency array *****

res = 12;                        % octave resolution
f_low = 30;                     % lower frequency bound
f_hi = 120;                     % upper frequency bound
n_octaves = log2(f_hi/f_low);   % number of octaves
f_bins = n_octaves * res + 1;   % number of frequency bins
f_bins_noalias = 19;            % frequency bins to avoid spatial aliasing.
% no aliasing
f_bins = f_bins_noalias;
f = zeros(1, f_bins);
f(1) = f_low;
% frequcny bins calculation.
for i = 2:f_bins                
   f(i) = f(i-1) * pow2(1/res);
end

% ***** monochromatic tensors initialisation *****
% init of the monochromatic pressure matrix
main_monochromatic_pressure = zeros(numel(y_control), ...
    numel(x_control), f_bins);
ANC_monochromatic_pressure = zeros(numel(y_control), ...
    numel(x_control), f_bins);
ANC_monochromatic_pressure_reg = zeros(numel(y_control), ...
    numel(x_control), f_bins);
plane_wave = zeros(numel(y_control), ...
    numel(x_control), f_bins);
% init of the 4D tensor
ANC_FFT_Data = zeros(numel(y_control), numel(x_control), f_bins, 3);

% ***** MAIN ARRAY SOUND FIELD *****
% propagation law is free-field monopole pressure p = A*e^(-jkr)/r

j = sqrt(-1);
c = 340;                                    % speed of sound

% ***** Per-frequency iteration *****
for i = 1:f_bins

    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    % ***** Matrix Calculation *****
    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_main)             % x-axis main array
                for q = 1:numel(y_main)         % y-axis main array
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                     r = sqrt((x_control(l) - x_main(p))^2 + ...
                         (y_control(m) - y_main(q))^2);
                 % pressure calculation as a summation of terms
                 % complex coefficient A is given by the (q,p)th source
                 % being on or off 
                 main_monochromatic_pressure(m,l,i) =  ... 
                     main_monochromatic_pressure(m,l,i) + ... 
                     main_array(q,p)*exp(-j*k*r)/r;
                end
            end
        end
    end
end

% ***** Plane Wave *****
% model reflections pl_wv = A*exp(-j[k_x*x+k_y*y])
% where k_x = kcos(theta), k_y = ksin(theta).

% angle = -30° with respect to the x axis.
theta_pw = -pi/6;
for i = 1:f_bins
    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber
    k_x(i) = k*cos(theta_pw);
    k_y(i) = k*sin(theta_pw);
    % amplitude coeff. (0.05)
    A = 0.05;
    for q = 1:numel(y_control)
        for p = 1:numel(x_control)
            plane_wave(q,p,i) = A*exp(-j*(k_x(i)*x_control(p) + ...
                k_y(i)*y_control(q))); 
        end
    end
end

% ***** Add reflections to the main sound field data *****
main_monochromatic_pressure = main_monochromatic_pressure + plane_wave;

% ***** NORMALIZATION *****
% get max value from the target matrix
% everytime shrinks dimensions
main_monochromatic_pressure_max = max(main_monochromatic_pressure);
main_monochromatic_pressure_max = max(main_monochromatic_pressure_max);
main_monochromatic_pressure_max = max(main_monochromatic_pressure_max);
main_monochromatic_pressure_max = abs(main_monochromatic_pressure_max);
% get min value from the target matrix
main_monochromatic_pressure_min = min(main_monochromatic_pressure);
main_monochromatic_pressure_min = min(main_monochromatic_pressure_min);
main_monochromatic_pressure_min = min(main_monochromatic_pressure_min);
main_monochromatic_pressure_min = abs(main_monochromatic_pressure_min);

% matrix normalisation
norm_main_monochromatic_pressure = (main_monochromatic_pressure - ...
    main_monochromatic_pressure_min)/(main_monochromatic_pressure_max - ...
    main_monochromatic_pressure_min);

% % normalisation for each frequency 
% norm_main_monochromatic_pressure = normalize(main_monochromatic_pressure, 3);

 for i = i:f_bins
    % mag
    ANC_FFT_Data(:,:,i,1) = abs(norm_main_monochromatic_pressure(:,:,i));
    % phase
    ANC_FFT_Data(:,:,i,2) = angle(norm_main_monochromatic_pressure(:,:,i));
    % 3rd empty channel
    ANC_FFT_Data(:,:,i,3) = zeros(numel(y_control), numel(x_control), 1);
 end
 
% ***** CNN PLACEMENT OPTIMIZATION *****
ANC_optimized_locations = predict(net, ANC_FFT_Data);

% ***** Locations hard clipping *****
ANC_optimized_locations_clipped = zeros(numel(ANC_optimized_locations), 1);
for i = 1:numel(ANC_optimized_locations)
    if ANC_optimized_locations(i) < 0.1
        ANC_optimized_locations_clipped(i) = 0.0;
    elseif ANC_optimized_locations(i) >= 0.1
        ANC_optimized_locations_clipped(i) = 1;
    end
end

ANC_locations_matrix = reshape(ANC_optimized_locations_clipped, ...
    numel(y_sources), numel(x_sources));

% ***** Line Disposition *****
secondary_array_reg = zeros(numel(y_sources), numel(x_sources));
for i = 1:numel(x_sources)
    for jj = 1:numel(y_sources)
        % sources on the first line -> ii == 1
        if i == 2 && (jj == 2 || jj == 4 || jj == 6 || jj == 8 || ...
                jj == 9 || jj == 11 || jj == 13 || jj == 15)
            secondary_array_reg(jj,i) = 1;
        else
            secondary_array_reg(jj,i) = 0;
        end
    end
end

% % ***** Random Disposition *****
% % Generate a symmetric disposition: halve the size of S, and all
% % subsequent random_array_half, rand_perm_array...
% % symmetry has to happen vertically on the x-axis, thus only y_sources
% % has to be halved.
% S = numel(x_sources)*numel(y_sources)/2;                       
% random_array_halved = zeros(S,1);       % initialize
% Kn = 4;
% % spread that many around
% % generate an array containing a random permutation of S elements in groups
% % of K. Size of the array is 1xK.
% rand_perm_array_halved = randperm(S,Kn);
% % elements of random_perm_array corresponding to the values of 
% % elements in random_array are given value 1. 
% rand_perm_array_halved = rand_perm_array_halved.';
% random_array_halved(rand_perm_array_halved) = ones(1,Kn);   
% % reshape from 1xS to size_y_sources/2 x size_x_sources matrix.
% secondary_array_reg = reshape(...
%     random_array_halved,[numel(y_sources)/2, numel(x_sources)]);
% % generate symmetrical matrix
% for i = 1:numel(y_sources)/2
%     j = numel(y_sources)/2 + i;
%     secondary_array_reg(j,:) = secondary_array_reg(...
%         numel(y_sources)/2 + 1 - i,:);
% end  

% reshaped option
secondary_array_reg1d = reshape(secondary_array_reg, ...
numel(y_sources)*numel(x_sources), 1);

% ***** PRESSURE MATCHING *****

% ***** Main Array Transfer Function *****
main_TF = zeros(numel(y_control), numel(x_control), f_bins);
for i = 1:f_bins
    
    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_main)             % x-axis main array
                for q = 1:numel(y_main)         % y-axis main array
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                     r = sqrt((x_control(l) - x_main(p))^2 + ...
                         (y_control(m) - y_main(q))^2);
                 % pressure calculation as a summation of terms
                 % complex coefficient A is given by the (q,p)th source
                 % being on or off 
                     main_TF(m,l,i) = main_TF(m,l,i) + ... 
                         main_array(q,p)*exp(-j*k*r)/r;
                end
            end
        end
    end
end
% reshape into a matrix -> each column corresponds to the TF at one freq.
main_TF = reshape(main_TF, numel(y_control)*numel(x_control), f_bins);

% ***** Secondary Array Transfer Function *****
% number of active speakers
K = 8;
secondary_TF_sparse = zeros(numel(y_control), numel(x_control), ...
    numel(y_sources), numel(x_sources), f_bins);
secondary_TF_sparse_reg = zeros(numel(y_control), numel(x_control), ...
    numel(y_sources), numel(x_sources), f_bins);
for i = 1:f_bins
    
    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_sources)          % x-axis secondary array
                for q = 1:numel(y_sources)      % y-axis secondary array
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                     r = sqrt((x_control(l) - x_sources(p))^2 + ...
                         (y_control(m) - y_sources(q))^2);
                 % pressure calculation as a summation of terms
                 % complex coefficient A is given by the (q,p)th source
                 % being on or off 
                     secondary_TF_sparse(m,l,q,p,i) =  ...
                         ANC_locations_matrix(q,p)*exp(-j*k*r)/r;
                     secondary_TF_sparse_reg(m,l,q,p,i) =  ...
                         secondary_array_reg(q,p)*exp(-j*k*r)/r;
                end
            end
        end
    end
end
secondary_TF_sparse = reshape(secondary_TF_sparse, ...
    numel(y_control)*numel(x_control), numel(y_sources)*numel(x_sources), f_bins);
secondary_TF_sparse_reg = reshape(secondary_TF_sparse_reg, ... 
    numel(y_control)*numel(x_control), numel(y_sources)*numel(x_sources), f_bins);

% calculate non-sparse secondary TF tensor 
secondary_TF = zeros(numel(y_control)*numel(x_control), K, f_bins);
secondary_TF_reg = zeros(numel(y_control)*numel(x_control), K, f_bins);
for i = 1:f_bins
    KK = 1;
    JJ = 1;
    for l = 1:numel(y_sources)*numel(x_sources)
        if ANC_optimized_locations_clipped(l) ~= 0
           secondary_TF(:,KK,i) = secondary_TF_sparse(:,l,i);
           KK = KK + 1;
        end
        if secondary_array_reg1d(l) ~= 0
           secondary_TF_reg(:,JJ,i) = secondary_TF_sparse_reg(:,l,i);
           JJ = JJ + 1;
        end
    end
end

% ***** Filters Calculation *****
filters = zeros(K, f_bins);
filters_reg = zeros(K, f_bins);
% regularisation parameter
delta = 1*10^-3;
% might change to match AE
delta_reg = 1.0*10^-3;
for i = 1:f_bins
    filters(:,i) = -inv(secondary_TF(:,:,i)'*secondary_TF(:,:,i) + ...
        delta*eye(K))*secondary_TF(:,:,i)'*main_TF(:,i);
    filters_reg(:,i) = -inv(secondary_TF_reg(:,:,i)'*secondary_TF_reg(:,:,i) + ...
        delta_reg*eye(K))*secondary_TF_reg(:,:,i)'*main_TF(:,i);
end
% transform to sparse matrix corresponding to all possible source locations
filters_sparse = zeros(numel(y_sources)*numel(x_sources), f_bins);
filters_sparse_reg = zeros(numel(y_sources)*numel(x_sources), f_bins);
for i = 1:f_bins
    KK = 1;
    JJ = 1;
    for l = 1:numel(x_sources)*numel(y_sources)
        if ANC_optimized_locations_clipped(l) ~= 0
            filters_sparse(l,i) = filters(KK,i);
            KK = KK + 1;
        end
        if secondary_array_reg1d(l) ~= 0
           filters_sparse_reg(l,i) = filters_reg(JJ,i);
           JJ = JJ + 1;
        end
    end
end
% reshape 
filters_sparse = reshape(filters_sparse, numel(y_sources), ...
    numel(x_sources), f_bins);
filters_sparse_reg = reshape(filters_sparse_reg, numel(y_sources), ...
    numel(x_sources), f_bins);

% ***** MAIN + SECONDARY ARRAY SOUND FIELD *****

% secondary array sound field
for i = 1:f_bins

    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    % ***** Matrix Calculation *****
    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_sources)          % x-axis secondary array
                for q = 1:numel(y_sources)      % y-axis secondary array
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                     r = sqrt((x_control(l) - x_sources(p))^2 + ...
                         (y_control(m) - y_sources(q))^2);
                 % pressure calculation as a summation of terms
                 % complex coefficient A is given by the (q,p)th source
                 % being on or off 
                 ANC_monochromatic_pressure(m,l,i) =  ... 
                     ANC_monochromatic_pressure(m,l,i) + ... 
                     filters_sparse(q,p,i)*exp(-j*k*r)/r;
                 ANC_monochromatic_pressure_reg(m,l,i) =  ... 
                     ANC_monochromatic_pressure_reg(m,l,i) + ... 
                     filters_sparse_reg(q,p,i)*exp(-j*k*r)/r;
                end
            end
        end
    end
end
% calculate rusulting sound fields
ANC_monochromatic_pressure = ANC_monochromatic_pressure + ...
    main_monochromatic_pressure;
ANC_monochromatic_pressure_reg = ANC_monochromatic_pressure_reg + ...
    main_monochromatic_pressure;

% calculate power
% calculate sound power (mean over freq) 
main_power = zeros(numel(y_control), numel(x_control));
ANC_power = zeros(numel(y_control), numel(x_control));
ANC_power_reg = zeros(numel(y_control), numel(x_control));
for l = 1:numel(x_control)
    for m = 1:numel(y_control)
        temp_sf = 0.0;
        temp_sf1 = 0.0;
        temp_sf2 = 0.0;
        for i = 1:f_bins
            temp_sf = temp_sf + ...
                norm(main_monochromatic_pressure(m,l,i))^2;
            temp_sf1 = temp_sf1 + ...
                norm(ANC_monochromatic_pressure(m,l,i))^2;
            temp_sf2 = temp_sf2 + ...
                norm(ANC_monochromatic_pressure_reg(m,l,i))^2;
        end
        main_power(m,l) = 10*log10(temp_sf/f_bins^2);
        ANC_power(m,l) = 10*log10(temp_sf1/f_bins^2);
        ANC_power_reg(m,l) = 10*log10(temp_sf2/f_bins^2);
    end
end
% power over frequency
main_power_f = zeros(1, f_bins);
ANC_power_f = zeros(1, f_bins);
ANC_power_f_reg = zeros(1, f_bins);
for i = 1:f_bins
    temp_sf = 0.0;
    temp_sf1 = 0.0;
    temp_sf2 = 0.0;
    for l = 1:numel(x_control)
        for m = 1:numel(y_control)
            temp_sf = temp_sf + norm(main_monochromatic_pressure(m,l,i))^2;
            temp_sf1 = temp_sf1 + norm(ANC_monochromatic_pressure(m,l,i))^2;
            temp_sf2 = temp_sf2 + norm(ANC_monochromatic_pressure_reg(m,l,i))^2;
        end
    end
    main_power_f(i) = temp_sf/(numel(x_control)*numel(y_control));
    ANC_power_f(i) = temp_sf1/(numel(x_control)*numel(y_control));
    ANC_power_f_reg(i) = temp_sf2/(numel(x_control)*numel(y_control));
end
% power average over freq
main_power_f_ave = 10*log10(sum(main_power_f)/f_bins);
ANC_power_f_ave = 10*log10(sum(ANC_power_f/f_bins));
ANC_power_f_ave_reg = 10*log10(sum(ANC_power_f_reg/f_bins));
% total power in db
main_power_f = 10*log10(main_power_f);
ANC_power_f = 10*log10(ANC_power_f);
ANC_power_f_reg = 10*log10(ANC_power_f_reg);
% save sound field power data
save(['main power ', ...
    num2str(x_main), '.mat'], 'main_power_f');
save(['ANC power ', ...
    num2str(x_main), '.mat'], 'ANC_power_f');
save(['ANC_reg power ', ...
    num2str(x_main), '.mat'], 'ANC_power_f_reg');

% power difference
ANC_power_diff = main_power - ANC_power;
ANC_power_diff_reg = main_power - ANC_power_reg;

% ***** METRICS *****
% ***** Insertion Loss *****
% IL = 10log10(||main_TF||^2/(||secondary_TF*filters + main_TF||)^2)
% Calculate average over freqs.
% reshape filters
filters_sparse = reshape(filters_sparse, numel(y_sources)*numel(x_sources), f_bins);
filters_sparse_reg = reshape(filters_sparse_reg, numel(y_sources)*numel(x_sources), f_bins);
IL = zeros(1, f_bins);
IL_reg = zeros(1, f_bins);
for i = 1:f_bins
    IL(i) = norm(main_TF(:,i))^2/ ... 
        norm(secondary_TF_sparse(:,:,i)*filters_sparse(:,i) + ...
        main_TF(:,i))^2;
    IL_reg(i) = norm(main_TF(:,i))^2/ ... 
        norm(secondary_TF_sparse_reg(:,:,i)*filters_sparse_reg(:,i) + ...
        main_TF(:,i))^2;
end
IL_ave = 10*log10(sum(IL)/f_bins);
IL_reg_ave = 10*log10(sum(IL_reg)/f_bins);
IL = 10*log10(IL);
IL_reg = 10*log10(IL_reg);
% save IL data
save(['IL ', ...
    num2str(x_main), '.mat'], 'IL');
save(['IL_reg ', ...
    num2str(x_main), '.mat'], 'IL_reg');

% ***** Array Effort *****
% AE = 10log10(norm(filters)^2/K)
% Calculate average over freq
AE = zeros(1, f_bins);
AE_reg = zeros(1, f_bins);
for i = 1:f_bins
    AE(i) = norm(filters(:,i))^2/K;
    AE_reg(i) = norm(filters_reg(:,i))^2/K;
end
AE_ave = 10*log10(sum(AE)/f_bins);
AE_reg_ave = 10*log10(sum(AE_reg)/f_bins);
AE = 10*log10(AE);
AE_reg = 10*log10(AE_reg);
% save AE data
save(['AE ', ...
    num2str(x_main), '.mat'], 'AE');
save(['AE_reg ', ...
    num2str(x_main), '.mat'], 'AE_reg');

% ***** PLOTS *****
% set alpha values
alphaVal = 0.6;
alphaVal_colorbar = 0.8;

% Plot average power levels
% add new row and column to avoid missing data when using surf
main_power(numel(y_control) + 1, :) = main_power(1,:);
main_power(:, numel(x_control) + 1) = main_power(:,1);
% main array columns and row are added to keep colormap low
ANC_power(numel(y_control) + 1, :) = main_power(1,1);
ANC_power(:, numel(x_control) + 1) = main_power(1,1);
ANC_power_reg(numel(y_control) + 1, :) = main_power(1,1);
ANC_power_reg(:, numel(x_control) + 1) = main_power(1,1);
ANC_power_diff(numel(y_control) + 1, :) = ANC_power_diff(1,:);
ANC_power_diff(:, numel(x_control) + 1) = ANC_power_diff(:,1);
ANC_power_diff_reg(numel(y_control) + 1, :) = ANC_power_diff_reg(1,:);
ANC_power_diff_reg(:, numel(x_control) + 1) = ANC_power_diff_reg(:,1);

% get max and min power values
max_power = max(main_power);
max_power = max(max_power);
min_power = min(main_power);
min_power = min(min_power) - 2;
min_max = [min_power max_power];
% get max and min ANC power values
% optimised array
max_power_ANC = max(ANC_power);
max_power_ANC = max(max_power_ANC);
min_power_ANC = min(ANC_power);
min_power_ANC = min(min_power_ANC);
% regular array
max_power_ANC_reg = max(ANC_power_reg);
max_power_ANC_reg = max(max_power_ANC_reg);
min_power_ANC_reg = min(ANC_power_reg);
min_power_ANC_reg = min(min_power_ANC_reg);
if max_power_ANC_reg > max_power_ANC
    max_power_ANC = max_power_ANC_reg;
end
if min_power_ANC_reg < min_power_ANC
    min_power_ANC = min_power_ANC_reg;
end

min_max_ANC = [min_power_ANC max_power_ANC];
%min_max_ANC = [-60 -30];

% min max power difference
% opt array
max_diff = max(ANC_power_diff);
max_diff = max(max_diff);
min_diff = min(ANC_power_diff);
min_diff = min(min_diff);
% reg array
max_diff_reg = max(ANC_power_diff_reg);
max_diff_reg = max(max_diff_reg);
min_diff_reg = min(ANC_power_diff_reg);
min_diff_reg = min(min_diff_reg);
if max_diff_reg > max_diff
    max_diff = max_diff_reg;
end
if min_diff_reg < min_diff
    min_diff_ANC = min_diff_reg;
end

min_max_diff = [min_diff max_diff];


% Optimised locations
% array containing values of all points on x-axis
x_source_points = zeros(numel(x_sources), 1);
for i=1:numel(x_sources)
    for jj=1:numel(y_sources)
      x_source_points((i-1)*numel(y_sources)+jj) = x_sources(i);  
    end
end
% array containing values of all points on y-axis
y_source_points = zeros(numel(y_sources), 1);
for i=1:numel(x_sources)
    for jj=1:numel(y_sources)
      y_source_points((i-1)*numel(y_sources)+jj) = ...
          y_sources(numel(y_sources)-jj+1);  
    end
end
% ''normalise'' locations array
% sound field power graph
ANC_optimized_locations_clipped_norm = zeros(numel(y_sources)*numel(x_sources), 1);
secondary_array_reg_norm = zeros(numel(y_sources)*numel(x_sources), 1);
main_array_norm = max_power*main_array;
for l = 1:numel(y_sources)*numel(x_sources)
        if ANC_optimized_locations_clipped(l) == 1
            ANC_optimized_locations_clipped_norm(l) = max_power_ANC;
        else
            ANC_optimized_locations_clipped_norm(l) = min_power_ANC;
        end
        if secondary_array_reg(l) == 1
            secondary_array_reg_norm(l) = max_power_ANC;
        else
            secondary_array_reg_norm(l) = min_power_ANC;
        end
        
end
% sound field differnce graph
ANC_optimized_locations_clipped_norm_diff = zeros(numel(y_sources)*numel(x_sources), 1);
secondary_array_reg_norm_diff = zeros(numel(y_sources)*numel(x_sources), 1);
main_array_norm_diff = 0.95*max_diff*main_array;
for l = 1:numel(y_sources)*numel(x_sources)
        if ANC_optimized_locations_clipped(l) == 1
            ANC_optimized_locations_clipped_norm_diff(l) = max_diff*0.95;
        else
            ANC_optimized_locations_clipped_norm_diff(l) = min_diff;
        end
        if secondary_array_reg(l) == 1
            secondary_array_reg_norm_diff(l) = max_diff*0.95;
        else
            secondary_array_reg_norm_diff(l) = min_diff;
        end   
end

figure(1)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
subplot(3,1,1)
% plot main array sound field
ss = surf(x_control_surf, y_control_surf, main_power,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Initial Sound Field Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max);
colormap(cool)
cc = colorbar;
view(0,90);
drawnow
% Make the colorbar transparent
cdata = cc.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal_colorbar * cdata(end,:));
cc.Face.Texture.ColorType = 'truecoloralpha';
cc.Face.Texture.CData = cdata;
drawnow
% Make sure that the renderer doesn't revert your changes
cc.Face.ColorBinding = 'discrete';
hold on
% main array
scatter3(x_scatter3,y_main,main_array_norm, 60, ...
    main_array_norm, 'filled', 'MarkerFaceAlpha', alphaVal);
hold off

% plot ANC sound fields 
% plot ANC sound field - optimised array
subplot(3,1,2)
ss = surf(x_control_surf, y_control_surf, ANC_power,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('ANC Sound Field Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('CNN Array', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_ANC);
colormap(cool)
cc = colorbar;
view(0,90);
drawnow
% Make the colorbar transparent
cdata = cc.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal_colorbar * cdata(end,:));
cc.Face.Texture.ColorType = 'truecoloralpha';
cc.Face.Texture.CData = cdata;
drawnow
% Make sure that the renderer doesn't revert your changes
cc.Face.ColorBinding = 'discrete';
hold on
% main array
scatter3(x_scatter3,y_main,main_array_norm, 60, ...
    main_array_norm, 'filled', 'MarkerFaceAlpha', alphaVal);
scatter3(x_source_points, y_source_points, ...
    ANC_optimized_locations_clipped_norm, 60, ...
    ANC_optimized_locations_clipped_norm, 'filled', 'MarkerFaceAlpha', ...
    alphaVal);
hold off

% plot ANC sound field - regular array
subplot(3,1,3)
ss = surf(x_control_surf, y_control_surf, ANC_power_reg,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('ANC Sound Field Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('Uniformly Spaced Array', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_ANC);
colormap(cool)
cc = colorbar;
view(0,90);
drawnow
% Make the colorbar transparent
cdata = cc.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal_colorbar * cdata(end,:));
cc.Face.Texture.ColorType = 'truecoloralpha';
cc.Face.Texture.CData = cdata;
drawnow
% Make sure that the renderer doesn't revert your changes
cc.Face.ColorBinding = 'discrete';
hold on
% main array
scatter3(x_scatter3,y_main,main_array_norm, 60, ...
    main_array_norm, 'filled', 'MarkerFaceAlpha', alphaVal);
scatter3(x_source_points, y_source_points, secondary_array_reg_norm, 60,...
    secondary_array_reg_norm, 'filled', 'MarkerFaceAlpha', alphaVal);
hold off

saveas(figure(1),...
    'Sound Field Power.png');

% plot sf power difference
figure(2)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
% plot power difference - optimised array
subplot(2,1,1)
ss = surf(x_control_surf, y_control_surf, ANC_power_diff,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('ANC Sound Field Power Difference [dB]', 'FontSize', 20, ...
    'FontWeight', 'normal', 'Interpreter', 'latex'); 
subtitle('CNN Array', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_diff);
colormap(copper)
cc = colorbar;
view(0,90);
drawnow
% Make the colorbar transparent
cdata = cc.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal_colorbar * cdata(end,:));
cc.Face.Texture.ColorType = 'truecoloralpha';
cc.Face.Texture.CData = cdata;
drawnow
% Make sure that the renderer doesn't revert your changes
cc.Face.ColorBinding = 'discrete';
hold on
% main array
scatter3(x_scatter3,y_main, main_array_norm_diff, 90, ...
    main_array_norm_diff, 'filled', 'MarkerFaceAlpha', 1);
scatter3(x_source_points, y_source_points, ...
    ANC_optimized_locations_clipped_norm_diff, 90, ...
    ANC_optimized_locations_clipped_norm_diff, 'filled', 'MarkerFaceAlpha', 1);
hold off

% plot ANC sound field - regular array
subplot(2,1,2)
ss = surf(x_control_surf, y_control_surf, ANC_power_diff_reg, ...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('ANC Sound Field Power Difference [dB]', 'FontSize', 20, ...
    'FontWeight', 'normal', 'Interpreter', 'latex'); 
subtitle('Uniformly Spaced Array', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_diff);
colormap(copper)
cc = colorbar;
view(0,90);
drawnow
% Make the colorbar transparent
cdata = cc.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal_colorbar * cdata(end,:));
cc.Face.Texture.ColorType = 'truecoloralpha';
cc.Face.Texture.CData = cdata;
drawnow
% Make sure that the renderer doesn't revert your changes
cc.Face.ColorBinding = 'discrete';
hold on
% main array
main_array_norm = max_power*main_array;
scatter3(x_scatter3,y_main, main_array_norm_diff, 90, ...
    main_array_norm_diff, 'filled', 'MarkerFaceAlpha', 1);
scatter3(x_source_points, y_source_points, secondary_array_reg_norm_diff, ...
    90, secondary_array_reg_norm_diff,'filled', 'MarkerFaceAlpha', 1);
hold off

saveas(figure(2),...
    'Sound Field Power Difference.png');

% plot optimised locations
figure(3)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
% plot nn optimized positions 
subplot(1,2,1)
scatter3(x_source_points, y_source_points, ANC_optimized_locations, 120, ...
    ANC_optimized_locations, 'filled', 'MarkerFaceAlpha', 0.8);
colormap(cool)
colorbar
axis equal
title('CNN Source Locations', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
view(0,90);                     % see the plot from above.
set(gcf,'color','w');           % set bg color to white.

subplot(1,2,2)
scatter3(x_source_points, y_source_points, ANC_optimized_locations_clipped, 120, ...
    ANC_optimized_locations_clipped, 'filled', 'MarkerFaceAlpha', 0.8);
colormap(cool)
colorbar
axis equal
title('CNN Source Locations [clipped]', 'FontSize', 20, ... 
    'FontWeight', 'normal', 'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
view(0,90);                     % see the plot from above.
set(gcf,'color','w');           % set bg color to white.
colorbar;

saveas(figure(3),...
    'Sources Placement.png');

% IL and AE
figure(4)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.

subplot(1,3,1)
plot(f, IL,'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
title('IL [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, IL_reg,'k*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
title('IL [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
hold off
legend(['CNN Array. Average IL = ' num2str(IL_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average IL = ' num2str(IL_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex');

subplot(1,3,2)
plot(f, AE,'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
title('AE [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, AE_reg,'k*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
hold off
legend(['CNN Array. Average AE = ' num2str(AE_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average AE = ' num2str(AE_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex');
subplot(1,3,3)
plot(f, main_power_f,'*-', 'Color', [0.20, 0.75, 0.90], ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 6, ...
    'linewidth', 1.5);
title('Sound Field Power', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Power Level [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, ANC_power_f,'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 6, 'linewidth', 1.5);
plot(f, ANC_power_f_reg,'k*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 6, 'linewidth', 1.5);
hold off
legend(['Initial SF Power. Average Power = ' num2str(main_power_f_ave) 'dB'], ...
    ['ANC CNN Array. Average Power = ' num2str(ANC_power_f_ave) 'dB'], ...
    ['ANC Uni. Array. Average Power = ' num2str(ANC_power_f_ave_reg) 'dB'], ...
    'location', 'NorthEast', 'FontSize', 11, 'Interpreter', 'latex');

saveas(figure(4),...
    'IL and AE.png');
