% ***** PLACEMENT OPTIMIZATION *****

% ***** Load Network *****
Select network from the Network folder
net = coder.loadDeepLearningNetwork( ...
    '');

% ***** Generate source and control region grids *****
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

% ***** monochromatic tensor initialisation *****
% init of the monochromatic pressure matrix
target_monochromatic_pressure = zeros(numel(y_control), ...
    numel(x_control), f_bins);

% ***** TARGET SOUND FIELD CALCULATION ***** 
% propagation law is free-field monopole pressure p = A*e^(-jkr)/r

j = sqrt(-1);
c = 340;                                    % speed of sound
% ***** Set active sources matrix for the desired sound field *****
% init
target_active_sources = zeros (numel(y_sources), numel(x_sources));

% % ***** Line Disposition *****
% for i = 1:numel(x_sources)
%     for jj = 1:numel(y_sources)
%         if i == 5 && (jj == 2 || jj == 4 || jj == 6 || jj == 8 || ...
%                 jj == 9 || jj == 11 || jj == 13 || jj == 15)
%             target_active_sources(jj,i) = 1;
%         else
%             target_active_sources(jj,i) = 0;
%         end
%     end
% end
% 
% % ***** Random Disposition *****
% % define symmetry
% symmetry_target = 1;
% % define number of active sources (1-valued elements in the grid)
% K = 8;
% if symmetry_target == 1
% K = K/2;
%  
% end   
% if symmetry_target == 0        % no symmetry
%     % S is a 1x2 array, make it scalar using the second element.                      
%     random_array = zeros(secondary_number,1);       % initialize  
%     % spread that many around
%     % generate an array containing a random permutation of S elements in groups
%     % of K. Size of the array is 1xK.
%     rand_perm_array = randperm(secondary_number,K);
%     % elements of random_perm_array corresponding to the values of 
%     % elements in random_array are given value 1. 
%     rand_perm_array = rand_perm_array.';
%     random_array(rand_perm_array) = ones(1,K);   
%     % reshape from 1xS to size_y_sources/2 x size_x_sources matrix.
%     target_active_sources = reshape(...
%         random_array,[numel(y_sources), numel(x_sources)]);
% elseif symmetry_target == 1    % y-axis symmetry
%     % Generate a symmetric disposition: halve the size of S, and all
%     % subsequent random_array_half, rand_perm_array...
%     % symmetry has to happen vertically on the x-axis, thus only y_sources
%     % has to be halved.
%     secondary_number = secondary_number/2;                       
%     random_array_halved = zeros(secondary_number,1);       % initialize
% 
%     % spread that many around
%     % generate an array containing a random permutation of S elements in groups
%     % of K. Size of the array is 1xK.
%     rand_perm_array_halved = randperm(secondary_number,K);
%     % elements of random_perm_array corresponding to the values of 
%     % elements in random_array are given value 1. 
%     rand_perm_array_halved = rand_perm_array_halved.';
%     random_array_halved(rand_perm_array_halved) = ones(1,K);   
%     % reshape from 1xS to size_y_sources/2 x size_x_sources matrix.
%     rand_perm_matrix_halved = reshape(...
%         random_array_halved,[numel(y_sources)/2, numel(x_sources)]);
%     % generate symmetrical matrix
%     target_active_sources = rand_perm_matrix_halved;
%     for i = 1:numel(y_sources)/2
%         jj = numel(y_sources)/2 + i;
%         target_active_sources(jj,:) = target_active_sources(...
%             numel(y_sources)/2 + 1 - i,:);
%     end
% end
% 
% % save distribution
% save('target_active_sources.mat',...
%     'target_active_sources'); 

% ***** load distribution *****
% loudspeaker dustributions are in the Target Distribution folder
n_fig = 1;
target_active_sources = load (['Target Distribution/target_sym' num2str(n_fig) '.mat']);
target_active_sources = struct2cell(target_active_sources);
target_active_sources = target_active_sources{1};
target_active_sources = reshape(...
    target_active_sources,[numel(y_sources), numel(x_sources)]);

% ***** Per-frequency iteration *****
for i = 1:f_bins

    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    % ***** Matrix Calculation *****
    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_sources)          % x-axis sources region
                for q = 1:numel(y_sources)      % y-axis sources region
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                 r = sqrt((x_control(l) - x_sources(p))^2 + ...
                     (y_control(m) - y_sources(q))^2);
                 % pressure calculation as a summation of terms
                 % being on or off 
                 % complex coefficient A is now 1.0.
                 target_monochromatic_pressure(m,l,i) =  ... 
                     target_monochromatic_pressure(m,l,i) + ... 
                     target_active_sources(q,p)*exp(-1j*k*r)/r;
                end
            end
        end
    end
end

% % ***** NORMALIZATION (needed if the network is trained on a normalised
% % dataset) *****
% 
% % get max value from the target matrix
% % everytime shrinks dimensions
% target_max = max(target_monochromatic_pressure);
% target_max = max(target_max);
% target_max = max(target_max);
% target_max = abs(target_max);
% % get min value from the target matrix
% target_min = min(target_monochromatic_pressure);
% target_min = min(target_min);
% target_min = min(target_min);
% target_min = abs(target_min);
% 
% % normalisation
% target_monochromatic_pressure = (target_monochromatic_pressure - ...
%     target_min)/(target_max - target_min);

% ***** Transform Raw_Data into FFT_Data having separate mag and phase.
% additional empty channel needed to treat FFT_Data as image data.
FFT_Data_target = zeros(numel(y_control), numel(x_control), f_bins, 3);

for i = 1:f_bins
    % mag
    % empirically scaled to get meaningful image data
    FFT_Data_target(:,:,i,1) = abs(20*log10(abs(...
        target_monochromatic_pressure(:,:,i))))/40;
    % phase
    % normalised to [0 1]
    FFT_Data_target(:,:,i,2) = angle(target_monochromatic_pressure(:,:,i))/pi + 0.5;
    % 3rd empty channel
    FFT_Data_target(:,:,i,3) = zeros(numel(y_control), numel(x_control), 1);
end

% ***** CNN PLACEMENT OPTIMIZATION *****
Optimized_locations = predict(net, FFT_Data_target);

% ***** Locations hard clipping *****
Optimized_locations_clipped = zeros(numel(Optimized_locations), 1);
for i = 1:numel(Optimized_locations)
    if Optimized_locations(i) < 0.32
        Optimized_locations_clipped(i) = 0.0;
    elseif Optimized_locations(i) >= 0.32
        Optimized_locations_clipped(i) = 1;
    end
end

% plot
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

% ***** OPTIMIZED SOUND FIELD CALCULATION *****
% resahpe optimized locations array
optimized_locations_matrix = reshape(Optimized_locations_clipped, ...
    numel(y_sources), numel(x_sources));

% init of the optmized monochromatic pressure matrix
optimized_monochromatic_pressure = zeros(numel(y_control), ...
    numel(x_control), f_bins);


% ***** Per-frequency iteration *****
% FIX this once aliasing problem is solved (array f is too long for plot)
f_bins = f_bins_noalias;
f_noalias = zeros(f_bins, 1);
for i = 1:f_bins
    f_noalias(i) = f(i);
end

for i = 1:f_bins

    lambda = c/f(i);                            
    k = 2*pi/lambda;                            % wavenumber

    % ***** Matrix Calculation *****
    for l = 1:numel(x_control)                  % x-axis control region
        for m = 1:numel(y_control)              % y-axis control region
            for p = 1:numel(x_sources)          % x-axis sources region
                for q = 1:numel(y_sources)      % y-axis sources region
                 % calculation of the element related to the i-th frequency bin 
                 % sound field given by p,q as the source position, control
                 % position given by l,m.
                 % distance calculation
                 r = sqrt((x_control(l) - x_sources(p))^2 + ...
                     (y_control(m) - y_sources(q))^2);
                 % pressure calculation as a summation of terms
                 % being on or off 
                 % complex coefficient A is now 1.0.
                 optimized_monochromatic_pressure(m,l,i) =  ... 
                     optimized_monochromatic_pressure(m,l,i) + ... 
                     optimized_locations_matrix(q,p)*exp(-1j*(k*r))/r;
                end
            end
        end
    end
end

% ***** NORMALIZATION *****
% get max value from the target matrix
% everytime shrinks dimensions
target_max = max(target_monochromatic_pressure);
target_max = max(target_max);
target_max = max(target_max);
target_max = abs(target_max);
% get min value from the target matrix
target_min = min(target_monochromatic_pressure);
target_min = min(target_min);
target_min = min(target_min);
target_min = abs(target_min);
% get max value from optimized matrix
optimized_max = max(optimized_monochromatic_pressure);
optimized_max = max(optimized_max);
optimized_max = max(optimized_max);
optimized_max = abs(optimized_max);
% get min value from optimized matrix
optimized_min = min(optimized_monochromatic_pressure);
optimized_min = min(optimized_min);
optimized_min = min(optimized_min);
optimized_min = abs(optimized_min);

% normalize matrices
target_monochromatic_pressure = (target_monochromatic_pressure - ...
    target_min)/(target_max - target_min);
optimized_monochromatic_pressure = (optimized_monochromatic_pressure - ...
    optimized_min)/(optimized_max - optimized_min);

% ***** ERROR METRICS CALCULATION *****

% SDR(pos) = 10log10(sum_f(|target_sf(pos,f)|)^2/sum_f(|target_sf(pos,f)-opt_sf(pos,f)|^2))
% MSE(pos) = 10log10(1/f_bins*sum_f(|target_sf(pos,f)-opt_sf(pos,f)|^2/|target_sf(pos,f)|^2))

% init error matrix
SDR_position = zeros(numel(y_control), numel(x_control));
MSE_position = zeros(numel(y_control), numel(x_control));
correlation_raw_pos = zeros(2, 2, numel(y_control), numel(x_control));
correlation_pos = zeros(numel(y_control), numel(x_control));

for l = 1:numel(x_control)
    for m = 1:numel(y_control)
        temp_error_num = 0.0;
        temp_error_den = 0.0;
        temp_ratio = 0.0;
        for i = 1:f_bins
            temp_error_num = temp_error_num + ...
                abs(target_monochromatic_pressure(m,l,i))^2;
            temp_error_den = temp_error_den + ...
                abs(target_monochromatic_pressure(m,l,i) - ...
                optimized_monochromatic_pressure(m,l,i))^2;
            temp_ratio = temp_ratio + (abs(target_monochromatic_pressure(m,l,i)...
                - optimized_monochromatic_pressure(m,l,i))^2/...
                abs(target_monochromatic_pressure(m,l,i))^2);
        end
        % SDR and MSE
        SDR_position(m,l) = 10*log10(temp_error_num/temp_error_den);
        MSE_position(m,l) = 10*log10(temp_ratio/f_bins); 
        % correlation_pos
        correlation_raw_pos(:, :, m, l) = ...
            corrcoef(target_monochromatic_pressure(m,l,:), ...
            optimized_monochromatic_pressure(m,l,:));
        correlation_pos(m,l) = correlation_raw_pos(1, 2, m, l);
    end
end

% save MSE_position
save(['MSE_position' num2str(n_fig) '.mat'], ...
    'MSE_position');

% SDR(f) = 10log10(sum_pos(|target_sf(pos,f)|)^2/sum_pos(|target_sf(pos,f)-opt_sf(pos,f)\^2))
% MSE(f) = 10log10(1/num_pos*sum_pos(|target_sf(pos,f)-opt_sf(pos,f)|^2/|target_sf(pos,f)|^2))

% init error matrix
SDR_frequency = zeros(f_bins, 1);
MSE_frequency = zeros(f_bins, 1);
MSE_frequency_wrong = zeros(f_bins, 1);

for i = 1:f_bins
    temp_error_num = 0.0;
    temp_error_den = 0.0;
    temp_ratio = 0.0;
    for m = 1:numel(y_control)
        for l = 1:numel(x_control)
            temp_error_num = temp_error_num + ...
                abs(target_monochromatic_pressure(m,l,i))^2;
            temp_error_den = temp_error_den + ...
                abs(target_monochromatic_pressure(m,l,i) - ...
                optimized_monochromatic_pressure(m,l,i))^2;
            temp_ratio = temp_ratio + (abs(target_monochromatic_pressure(m,l,i)...
                - optimized_monochromatic_pressure(m,l,i))^2/...
                abs(target_monochromatic_pressure(m,l,i))^2);
        end
    end
    SDR_frequency(i) = 10*log10(temp_error_num/temp_error_den);
    MSE_frequency_wrong(i) = 10*log10(temp_error_den/(temp_error_num*...
        numel(x_control)*numel(y_control)));
    MSE_frequency(i) = 10*log10(temp_ratio/(numel(x_sources)*numel(y_sources)));
end

% save MSE_frequency
save(['MSE_frequency' num2str(n_fig) '.mat'], ...
    'MSE_frequency');

% correlation_f
correlation_raw_f = zeros(2, 2, f_bins);
correlation_f = zeros(1, f_bins);
for i = 1:f_bins
    correlation_raw_f(:,:,i) = corrcoef(target_monochromatic_pressure(:,:,i), ...
        optimized_monochromatic_pressure(:,:,i));
    correlation_f(i) = correlation_raw_f(1,2,i);
end

% save correlation data
save(['Corr_freq_' num2str(n_fig) '.mat'], ...
    'correlation_f');
save(['Corr_pos_' num2str(n_fig) '.mat'], ...
    'correlation_pos');

% ***** PLOTS *****

% set alpha values
alphaVal = 0.6;
alphaVal_colorbar = 0.8;

figure(1)
% make figure as large as screen size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% plot active sources for the target field
subplot(1,3,1)
target_active_sources = reshape(target_active_sources, ...
    [numel(y_sources)*numel(x_sources) 1]);
scatter3(x_source_points, y_source_points, target_active_sources, 120, ...
    target_active_sources, 'filled', 'MarkerFaceAlpha', 0.8);
colormap(cool)
colorbar;
axis equal;
view(0,90);
title('Active Target Sources', ...
     'FontSize', 20, 'FontWeight', 'normal', 'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);

% plot nn optimized positions 
subplot(1,3,2)
scatter3(x_source_points, y_source_points, Optimized_locations, 120, ...
    Optimized_locations, 'filled', 'MarkerFaceAlpha', 0.8);
colormap(cool)
% flip colormap
% oldcmap = colormap;
% colormap(flipud(oldcmap));
colorbar
axis equal
title('CNN Source Locations', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
view(0,90);                     % see the plot from above.
set(gcf,'color','w');           % set bg color to white.

subplot(1,3,3)
scatter3(x_source_points, y_source_points, Optimized_locations_clipped, 120, ...
    Optimized_locations_clipped, 'filled', 'MarkerFaceAlpha', 0.8);
colormap(cool)
% flip colormap
% oldcmap = colormap;
% colormap(flipud(oldcmap));
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

saveas(figure(1), ...
    'Opt Sources.png');

figure(2)
% make figure as large as screen size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.

subplot(1,4,1)
% Plot average power levels
% calculate mean sound power (mean over freq) 
target_freq_mean_power = zeros(numel(y_control) + 1, numel(x_control) + 1);
for l = 1:numel(x_control)
    for m = 1:numel(y_control)
        temp_sf = 0.0;
        for i = 1:f_bins
            temp_sf = temp_sf + ...
                abs(target_monochromatic_pressure(m,l,i)^2);
        end
        target_freq_mean_power(m,l) = 10*log10(temp_sf/f_bins^2);
    end
end
% add new row and column to avoid missing data when using surf
target_freq_mean_power(numel(y_control) + 1, :) = ...
    target_freq_mean_power(1,:);
target_freq_mean_power(:, numel(x_control) + 1) = ...
    target_freq_mean_power(:,1);

ss = surf(x_control_surf, y_control_surf, target_freq_mean_power,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Target SF Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
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

subplot(1,4,2)
optimized_freq_mean_power = zeros(numel(y_control) + 1, numel(x_control) + 1);
for l = 1:numel(x_control)
    for m = 1:numel(y_control)
        temp_sf = 0.0;
        for i = 1:f_bins
           temp_sf = temp_sf + ...
                abs(optimized_monochromatic_pressure(m,l,i)^2);
        end
        optimized_freq_mean_power(m,l) = 10*log10(temp_sf/f_bins^2);
    end
end
% add new row and column to avoid missing data when using surf
optimized_freq_mean_power(numel(y_control) + 1, :) = ...
    optimized_freq_mean_power(1,:);
optimized_freq_mean_power(:, numel(x_control) + 1) = ...
    optimized_freq_mean_power(:,1);

ss = surf(x_control_surf, y_control_surf, optimized_freq_mean_power, ...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('CNN SF Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
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

% add new row and column to avoid missing data when using surf
SDR_position_surf = SDR_position;
SDR_position_surf(numel(y_control) + 1, :) = SDR_position_surf(1,:);
SDR_position_surf(:, numel(x_control) + 1) = SDR_position_surf(:,1);

MSE_position_surf = MSE_position;
MSE_position_surf(numel(y_control) + 1, :) = MSE_position_surf(1,:);
MSE_position_surf(:, numel(x_control) + 1) = MSE_position_surf(:,1);

correlation_pos_surf = correlation_pos;
correlation_pos_surf(numel(y_control) + 1, :) = correlation_pos_surf(1,:);
correlation_pos_surf(:, numel(x_control) + 1) = correlation_pos_surf(:,1);

% plot MSE_position
subplot(1,4,3)
ss = surf(x_control_surf, y_control_surf, MSE_position_surf,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('MSE [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
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

% plot correlation_pos
subplot(1,4,4)
ss = surf(x_control_surf, y_control_surf, real(correlation_pos_surf),...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
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
% save
saveas(figure(2), ...
    'corr pos.png');

% plot MSE_frequency

figure(3)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
subplot(1,2,1)
plot(f_noalias, MSE_frequency,'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
title('MSE(f) [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
% plot correlation_f
subplot(1,2,2)
plot(f_noalias, real(correlation_f),'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
title('Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('{Re\{.\}}', 'Interpreter', 'latex', 'FontSize', 15);
grid on
% subplot(1,3,3)
% plot(f_noalias, imag(correlation_f),'m*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
%     'MarkerSize', 8, 'linewidth', 1.5);
% title('Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
%     'Interpreter', 'latex'); 
% subtitle('', 'FontSize', 16)
% xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
% ylabel('{Im\{.\}}', 'Interpreter', 'latex', 'FontSize', 15);
% grid on
saveas(figure(3), ...
    'corr freq.png');

