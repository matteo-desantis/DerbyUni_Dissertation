% ***** DATA GENERATION *****

% ***** Generate grids for sources and control region *****
% Sources
x_sources = 0.0:1.0:5.0;
y_sources = -7.5:1.0:7.5;

% Control region
x_control = 20.0:2:30.0;
y_control = -10.0:2:10.0;

% ***** Main Array Grid *****
x_main = -50;
y_main = y_sources;
main_array = ones(numel(y_main), 1);
x_scatter3 = -50*ones(numel(y_main), 1);

% ***** 1/12 octave spaced frequency array *****

res = 12;                       % octave resolution
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

% ***** Set number of samples to be generated *****
n_samples = 5000;
noise_iterations = 5;


% ***** Y-axis SYMMETRY *****
% symmetry = 0 -> no symmetry
% symmetry = 1 -> y-axis sources symmetry 
symmetry = 0;
% define number of active sources (1-valued elements in the grid)
K = 8;
if symmetry == 1
    K = K/2;
end

% ***** ITERATE *****
for iteration = 1:n_samples
    

    % ***** Random permutations of k active sources each time *****

    % sample size
    size_x_sources = size(x_sources);
    size_y_sources = size(y_sources);
    S = size_x_sources.*size_y_sources;
    
    if symmetry == 0        % no symmetry
        % S is a 1x2 array, make it scalar using the second element. 
        S = S(2);                       
        random_array = zeros(S,1);       % initialize
        
        % spread that many around
        % generate an array containing a random permutation of S elements in groups
        % of K. Size of the array is 1xK.
        rand_perm_array = randperm(S,K);
        % elements of random_perm_array corresponding to the values of 
        % elements in random_array are given value 1. 
        rand_perm_array = rand_perm_array.';
        random_array(rand_perm_array) = ones(1,K);   
        % reshape from 1xS to size_y_sources/2 x size_x_sources matrix.
        rand_perm_matrix = reshape(...
            random_array,[size_y_sources(2), size_x_sources(2)]);
    
    elseif symmetry == 1    % y-axis symmetry
        % Generate a symmetric disposition: halve the size of S, and all
        % subsequent random_array_half, rand_perm_array...
        % symmetry has to happen vertically on the x-axis, thus only y_sources
        % has to be halved.
        S = S(2)/2;                       
        random_array_halved = zeros(S,1);       % initialize

        % spread that many around
        % generate an array containing a random permutation of S elements in groups
        % of K. Size of the array is 1xK.
        rand_perm_array_halved = randperm(S,K);
        % elements of random_perm_array corresponding to the values of 
        % elements in random_array are given value 1. 
        rand_perm_array_halved = rand_perm_array_halved.';
        random_array_halved(rand_perm_array_halved) = ones(1,K);   
        % reshape from 1xS to size_y_sources/2 x size_x_sources matrix.
        rand_perm_matrix_halved = reshape(...
            random_array_halved,[size_y_sources(2)/2, size_x_sources(2)]);
        % generate symmetrical matrix
        rand_perm_matrix = rand_perm_matrix_halved;
        for i = 1:size_y_sources(2)/2
            j = (size_y_sources(2))/2 + i;
            rand_perm_matrix(j,:) = rand_perm_matrix(...
                size_y_sources(2)/2 + 1 - i,:);
        end
        % flatten to get label array
        random_array = reshape(rand_perm_matrix, ...
            size_y_sources(2)*size_x_sources(2), 1);   
    end
    
    % ***** monochromatic tensor initialisation *****
        % init of the monochromatic pressure matrix
        monochromatic_pressure = zeros(numel(y_control), ...
            numel(x_control), f_bins);

%     % ***** Plot *****
%     % array containing values of all points on x-axis
%     for i = 1:size_x_sources(2)
%         for j = 1:size_y_sources(2)
%           x_source_points((i-1)*size_y_sources(2)+j) = x_sources(i);  
%         end
%     end
%     % array containing values of all points on y-axis
%     for i = 1:size_x_sources(2)
%         for j = 1:size_y_sources(2)
%           y_source_points((i-1)*size_y_sources(2)+j) = ...
%               y_sources(size_y_sources(2)-j+1);  
%         end
%     end
% 
%     figure(1)
%     % figure bigger
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     set(gcf,'color','w');           % set bg color to white.
%     % secondary array
%     scatter3(x_source_points,y_source_points,random_array.', 80, ...
%         random_array.', 'filled', 'MarkerFaceAlpha', 0.8);
%     hold on
%     % main array
%     scatter3(x_scatter3,y_main,main_array.', 80, ...
%         main_array.', 'filled', 'MarkerFaceAlpha', 0.8);
%     hold off
%     colormap(cool)
%     % flip colormap
%     % oldcmap = colormap;
%     % colormap(flipud(oldcmap));
%     colorbar
%     axis equal
%     title('Sources Locations', 'FontSize', 20, 'FontWeight', 'normal'); 
%     subtitle('', 'FontSize', 16)
%     xlabel('x [m]');
%     ylabel('y [m]');
%     view(0,90);                     % see the plot from above.
% %     % saving
% %     saveas(gcf, ...
% %     'Training Sources.png');
         
    % ***** SOUND FIELD CALCULATION ***** 
    % propagation law is free-field monopole pressure p = A*e^(-jkr)/r

    j = sqrt(-1);
    c = 340;                                    % speed of sound

    % init of the 4D tensor for storage
    FFT_Data = zeros(numel(y_control), numel(x_control), f_bins, 3);

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
                     % complex coefficient A is given by the (q,p)th source
                     % being on or off 
                     monochromatic_pressure(m,l,i) =  ... 
                         monochromatic_pressure(m,l,i) + ... 
                         rand_perm_matrix(q,p)*exp(-j*k*r)/r;
                    end
                end
            end
        end

%         % ***** Plotting *****
%         % change figure size
%         pos = get(gcf, 'Position');     % gives x left, y bottom, width, height
%         width = pos(3);
%         height = pos(4);
%         figure('Renderer', 'painters', 'Position', [pos(1) pos(2) 1.60*pos(3) 1.60*pos(4)]);
%         % mag
%         subplot(1,2,1)
%         s = surf(x_control,y_control,20*log10(abs(monochromatic_pressure(:,:,i))), ...
%             'FaceAlpha',0.5); %'FaceColor', 'interp');
%         s.EdgeColor = 'none';
%         axis equal;
%         title('Mag [dB]', 'FontSize', 20, 'FontWeight', 'normal'); 
%         subtitle(['Frequency = ' num2str(f(i)) ' [Hz]'], 'FontSize', 16, 'FontWeight', 'normal')
%         xlabel('x [m]');
%         ylabel('y [m]');
%         colormap(cool)
%         colorbar;
%         view(0,90);
%         % phase
%         subplot(1,2,2)
%         angle_rad = angle(monochromatic_pressure(:,:,i));
%         angle_deg = angle_rad/(2*pi)*180;
%         s = surf(x_control,y_control,angle_deg,'FaceAlpha',0.5); % 'FaceColor', 'interp');
%         s.EdgeColor = 'none';
%         axis equal;
%         title('Phase [°]', 'FontSize', 20, 'FontWeight', 'normal'); 
%         subtitle(['Frequency = ' num2str(f(i)) ' [Hz]'], 'FontSize', 16, 'FontWeight', 'normal')
%         xlabel('x [m]');
%         xlabel('x [m]');
%         ylabel('y [m]');
%         colormap(cool)
%         colorbar;
%         view(0,90);
%         % set bg color to white
%         set(gcf,'color','w');
%         % saving
%         current_image = sprintf(...
%         'ControlANC %f.png', f(i));
%         saveas(gcf, current_image);

%     % ***** ADD MAIN ARRAY *****
%     % ***** Per-frequency iteration *****
% 
%         % ***** Matrix Calculation *****
%         for l = 1:numel(x_control)                  % x-axis control region
%             for m = 1:numel(y_control)              % y-axis control region
%                 for p = 1:numel(x_main)             % x-axis main array
%                     for q = 1:numel(y_main)         % y-axis main array
%                      % calculation of the element related to the i-th frequency bin 
%                      % sound field given by p,q as the source position, control
%                      % position given by l,m.
%                      % distance calculation
%                          r = sqrt((x_control(l) - x_main(p))^2 + ...
%                              (y_control(m) - y_main(q))^2);
%                      % pressure calculation as a summation of terms
%                      % complex coefficient A is given by the (q,p)th source
%                      % being on or off 
%                      monochromatic_pressure(m,l,i) =  ... 
%                          monochromatic_pressure(m,l,i) + ... 
%                          main_array(q,p)*exp(-j*k*r)/r;
%                     end
%                 end
%             end
%         end
    
    end
    
    % ***** NORMALIZATION *****
    % get max value from the target matrix
    % everytime shrinks dimensions
    monochromatic_pressure_max = max(monochromatic_pressure);
    monochromatic_pressure_max = max(monochromatic_pressure_max);
    monochromatic_pressure_max = max(monochromatic_pressure_max);
    monochromatic_pressure_max = abs(monochromatic_pressure_max);
    % get min value from the target matrix
    monochromatic_pressure_min = min(monochromatic_pressure);
    monochromatic_pressure_min = min(monochromatic_pressure_min);
    monochromatic_pressure_min = min(monochromatic_pressure_min);
    monochromatic_pressure_min = abs(monochromatic_pressure_min);
    
    % overall matrix normalisation
    monochromatic_pressure = (monochromatic_pressure - ...
        monochromatic_pressure_min)/(monochromatic_pressure_max - ...
        monochromatic_pressure_min);

%     % normalisation for each frequency 
%     monochromatic_pressure = normalize(monochromatic_pressure, 3);

%     % overall zero mean 1 std
%     monochromatic_pressure = reshape(monochromatic_pressure, ...
%         numel(y_control)*numel(x_control)*f_bins, 1);
%     monochromatic_pressure = normalize(monochromatic_pressure);
%     monochromatic_pressure = reshape(monochromatic_pressure, ...
%         numel(y_control), numel(x_control), f_bins);
    
    % ***** Store sound field data as a 3d image (4D tensor) *****
    for i = i:f_bins
        % mag
        % empirically scaled to get meaningful data
        % FFT_Data(:,:,i,1) = abs(20*log10(abs(...
        %     monochromatic_pressure(:,:,i))))/40;
        FFT_Data(:,:,i,1) = abs(monochromatic_pressure(:,:,i));
        % phase
        % normalised to [0 1]
        % FFT_Data(:,:,i,2) = angle(monochromatic_pressure(:,:,i))/pi + 0.5;
        FFT_Data(:,:,i,2) = angle(monochromatic_pressure(:,:,i));
        % 3rd empty channel
        FFT_Data(:,:,i,3) = zeros(numel(y_control), numel(x_control), 1);
    end
    
    % ***** ADD NOISE *****
    % white gaussian noise added to the image data
    % each image is being noised noise_iterations times.
    for noise_iteration = 1:noise_iterations 
        % include uncorrupted data.
        if noise_iteration == 1
            FFT_Data_Noise = FFT_Data;
        else
        % third argument specifies variance. Average between noise and non
        % noise to reduce variance.
            FFT_Data_Noise = FFT_Data/2 + ...
                imnoise(FFT_Data, 'gaussian', 0.01)/2;
        end
        % ***** Saving files *****
        % joint example and label identifier
        current_counter = (iteration - 1)*noise_iterations ...
            + noise_iteration;

        % save image as a .mat file, while updating image data name
        save(['' ...
            num2str(current_counter) '.mat'], 'FFT_Data_Noise');
        
        % save labels as a .mat file
        save(['' ...
            num2str(current_counter) '.mat'], 'random_array');
        
%         figure(current_counter + 1)
%         % show generated image
%         current_image(:,:,:,1) = FFT_Data_Noise(:,:,1,:);
%         imshow(current_image, 'InitialMagnification',10000);


    end
 
end

