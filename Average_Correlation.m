% ***** LOAD DATA *****
% data is stored in Correlation and MSE Folders

% ***** CNN *****
% init
n_fields = 15;
cnn_correlation_f = zeros(19, n_fields);
cnn_correlation_pos = zeros(11, 6, n_fields);
cnn_MSE_f = zeros(19, n_fields);
cnn_MSE_pos = zeros(11, 6, n_fields);
for i = 1:n_fields
    cnn_correlation_f(:,i) = importdata(['Correlation/CNN_4D_024/Corr_freq_' num2str(i) '.mat']);
    cnn_correlation_pos(:,:,i) = importdata(['Correlation/CNN_4D_024/Corr_pos_' num2str(i) '.mat']);
    cnn_MSE_f(:,i) = importdata(['MSE/CNN_4D_024/MSE_frequency' num2str(i) '.mat']);
    cnn_MSE_pos(:,:,i) = importdata(['MSE/CNN_4D_024/MSE_position' num2str(i) '.mat']);
end

% ***** CNN2 Group 2 *****
cnn2_g2_correlation_f = zeros(19, n_fields);
cnn2_g2_correlation_pos = zeros(11, 6, n_fields);
cnn2_g2_MSE_f = zeros(19, n_fields);
cnn2_g2_MSE_pos = zeros(11, 6, n_fields);
for i = 1:n_fields
    cnn2_g2_correlation_f(:,i) = importdata(['Correlation/CNN_4D_G2_026/Corr_freq_' num2str(i) '.mat']);
    cnn2_g2_correlation_pos(:,:,i) = importdata(['Correlation/CNN_4D_G2_026/Corr_pos_' num2str(i) '.mat']);
    cnn2_g2_MSE_f(:,i) = importdata(['MSE/CNN_4D_G2_026/MSE_frequency' num2str(i) '.mat']);
    cnn2_g2_MSE_pos(:,:,i) = importdata(['MSE/CNN_4D_G2_026/MSE_position' num2str(i) '.mat']);
end

% ***** CNN group 2 *****
cnn_g2_correlation_f = zeros(19, n_fields);
cnn_g2_correlation_pos = zeros(11, 6, n_fields);
cnn_g2_MSE_f = zeros(19, n_fields);
cnn_g2_MSE_pos = zeros(11, 6, n_fields);
for i = 1:n_fields
    cnn_g2_correlation_f(:,i) = importdata(['Correlation/CNN_4D_G2/Corr_freq_' num2str(i) '.mat']);
    cnn_g2_correlation_pos(:,:,i) = importdata(['Correlation/CNN_4D_G2/Corr_pos_' num2str(i) '.mat']);
    cnn_g2_MSE_f(:,i) = importdata(['MSE/CNN_4D_G2/MSE_frequency' num2str(i) '.mat']);
    cnn_g2_MSE_pos(:,:,i) = importdata(['MSE/CNN_4D_G2/MSE_position' num2str(i) '.mat']);
end

% ***** CALCULATE METRICS *****
% ***** Mean *****
% cnn
cnn_correlation_f_ave = sum(cnn_correlation_f, 2)/n_fields;
cnn_correlation_f_glob_ave  = sum(cnn_correlation_f_ave)/f_bins;
cnn_correlation_pos_ave = sum(cnn_correlation_pos, 3)/n_fields;
% eliminate -inf data for the MSE calc
cnn_MSE_f_ave = 0;
cnn_MSE_pos_ave = 0;
for i = 1:n_fields
    if i ~= 5
        % undo log 
        cnn_MSE_f_ave = cnn_MSE_f_ave + 10.^(cnn_MSE_f(:,i)/10);
        cnn_MSE_pos_ave = cnn_MSE_pos_ave + 10.^(cnn_MSE_pos(:,:,i)/10);
    end
end
cnn_MSE_f_ave = cnn_MSE_f_ave/(n_fields - 1);
cnn_MSE_f_glob_ave = 10*log10(sum(cnn_MSE_f_ave)/f_bins);
cnn_MSE_f_ave = 10*log10(cnn_MSE_f_ave);
cnn_MSE_pos_ave = 10*log10(cnn_MSE_pos_ave/(n_fields - 1));

% cnn2 group 2
cnn2_g2_correlation_f_ave = sum(cnn2_g2_correlation_f, 2)/n_fields;
cnn2_g2_correlation_pos_ave = sum(cnn2_g2_correlation_pos, 3)/n_fields;
cnn2_g2_correlation_f_glob_ave  = sum(cnn2_g2_correlation_f_ave)/f_bins;
% eliminate -inf data for the MSE calc
cnn2_g2_MSE_f_ave = 0;
cnn2_g2_MSE_pos_ave = 0;
for i = 1:n_fields
    if i ~= 5
        % undo log
        cnn2_g2_MSE_f_ave = cnn2_g2_MSE_f_ave + 10.^(cnn2_g2_MSE_f(:,i)/10);
        cnn2_g2_MSE_pos_ave = cnn2_g2_MSE_pos_ave + ...
            10.^(cnn2_g2_MSE_pos(:,:,i)/10);
    end
end
cnn2_g2_MSE_f_ave = cnn2_g2_MSE_f_ave/(n_fields - 1);
cnn2_g2_MSE_f_glob_ave = 10*log10(sum(cnn2_g2_MSE_f_ave)/f_bins);
cnn2_g2_MSE_f_ave = 10*log10(cnn2_g2_MSE_f_ave);
cnn2_g2_MSE_pos_ave = 10*log10(cnn2_g2_MSE_pos_ave/(n_fields - 1));

% cnn group 2
cnn_g2_correlation_f_ave = sum(cnn_g2_correlation_f, 2)/n_fields;
cnn_g2_correlation_pos_ave = sum(cnn_g2_correlation_pos, 3)/n_fields;
cnn_g2_correlation_f_glob_ave  = sum(cnn_g2_correlation_f_ave)/f_bins;
% eliminate -inf data for the MSE calc
cnn_g2_MSE_f_ave = 0;
cnn_g2_MSE_pos_ave = 0;
for i = 1:n_fields
    if i ~=5
        % undo log
        cnn_g2_MSE_f_ave = cnn_g2_MSE_f_ave + 10.^(cnn_g2_MSE_f(:,i)/10);
        cnn_g2_MSE_pos_ave = cnn_g2_MSE_pos_ave + ...
            10.^(cnn_g2_MSE_pos(:,:,i)/10);
    end
end
cnn_g2_MSE_f_ave = cnn_g2_MSE_f_ave/(n_fields - 1);
cnn_g2_MSE_f_glob_ave = 10*log10(sum(cnn_g2_MSE_f_ave)/f_bins);
cnn_g2_MSE_f_ave = 10*log10(cnn_g2_MSE_f_ave);
cnn_g2_MSE_pos_ave = 10*log10(cnn_g2_MSE_pos_ave/(n_fields - 1));

% ***** standard deviation (correlation_f) *****
cnn_correlation_f_std = std(real(cnn_correlation_f), 0, 2);
cnn2_g2_correlation_f_std = std(real(cnn2_g2_correlation_f), 0, 2);
cnn_g2_correlation_f_std = std(real(cnn_g2_correlation_f), 0, 2);

% ***** median (correlation_f && MSE_f) *****
cnn_correlation_f_med = median(real(cnn_correlation_f), 2);
cnn2_g2_correlation_f_med = median(real(cnn2_g2_correlation_f), 2);
cnn_g2_correlation_f_med = median(real(cnn_g2_correlation_f), 2);
cnn_MSE_f_med = median(real(cnn_MSE_f), 2);
cnn2_g2_MSE_f_med = median(real(cnn2_g2_MSE_f), 2);
cnn_g2_MSE_f_med = median(real(cnn_g2_MSE_f), 2);

% ***** trimmed mean (correlation_f) *****
cnn_correlation_f_trm = trimmean(real(cnn_correlation_f), 10, 'round', 2);
cnn2_g2_correlation_f_trm = trimmean(real(cnn2_g2_correlation_f), 10, 'round', 2);
cnn_g2_correlation_f_trm = trimmean(real(cnn_g2_correlation_f), ...
    10, 'round', 2);

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

% control point grid for surf usage
x_control_surf = 19.0:2:31.0;
y_control_surf = -11.0:2:11.0;

% get max and min values of correlation_pos
% max
max_correlation = max(cnn_correlation_pos_ave);
max_correlation = real(max(max_correlation));
temp_max_correlation = max(cnn2_g2_correlation_pos_ave);
temp_max_correlation = real(max(temp_max_correlation));
if temp_max_correlation > max_correlation
    max_correlation = temp_max_correlation;
end
temp_max_correlation = max(cnn_g2_correlation_pos_ave);
temp_max_correlation = real(max(temp_max_correlation));
if temp_max_correlation > max_correlation
    max_correlation = temp_max_correlation;
end
% min
min_correlation = min(cnn_correlation_pos_ave);
min_correlation = real(min(min_correlation));
temp_min_correlation = min(cnn2_g2_correlation_pos_ave);
temp_min_correlation = real(min(temp_min_correlation));
if temp_min_correlation < min_correlation
    min_correlation = temp_min_correlation;
end
temp_min_correlation = min(cnn_g2_correlation_pos_ave);
temp_min_correlation = real(min(temp_min_correlation));
if temp_min_correlation < min_correlation
    min_correlation = temp_min_correlation;
end

min_max_corr = [min_correlation max_correlation];

% get max and min values of MSE_pos
% max
max_MSE = max(cnn_MSE_pos_ave);
max_MSE = real(max(max_MSE));
temp_max_MSE = max(cnn2_g2_MSE_pos_ave);
temp_max_MSE = real(max(temp_max_MSE));
if temp_max_MSE > max_MSE
    max_MSE = temp_max_MSE;
end
temp_max_MSE = max(cnn_g2_MSE_pos_ave);
temp_max_MSE = real(max(temp_max_MSE));
if temp_max_MSE > max_MSE
    max_MSE = temp_max_MSE;
end
% min
min_MSE = min(cnn_MSE_pos_ave);
min_MSE = real(min(min_MSE));
temp_min_MSE = min(cnn2_g2_MSE_pos_ave);
temp_min_MSE = real(min(temp_min_MSE));
if temp_min_MSE < min_MSE
    min_MSE = temp_min_MSE;
end
temp_min_MSE = min(cnn_g2_MSE_pos_ave);
temp_min_MSE = real(min(temp_min_MSE));
if temp_min_MSE < min_MSE
    min_MSE = temp_min_MSE;
end

min_max_MSE = [min_MSE max_MSE];

% ***** PLOT *****
% add new row and column to avoid missing data when using surf
cnn_correlation_pos_ave(numel(y_control) + 1, :) = ...
    cnn_correlation_pos_ave(1,:);
cnn_correlation_pos_ave(:, numel(x_control) + 1) = ...
    cnn_correlation_pos_ave(:,1);
cnn2_g2_correlation_pos_ave(numel(y_control) + 1, :) = ...
    cnn2_g2_correlation_pos_ave(1,:);
cnn2_g2_correlation_pos_ave(:, numel(x_control) + 1) = ...
    cnn2_g2_correlation_pos_ave(:,1);
cnn_g2_correlation_pos_ave(numel(y_control) + 1, :) = ...
    cnn_g2_correlation_pos_ave(1,:);
cnn_g2_correlation_pos_ave(:, numel(x_control) + 1) = ...
    cnn_g2_correlation_pos_ave(:,1);
cnn_MSE_pos_ave(numel(y_control) + 1, :) = ...
    cnn_MSE_pos_ave(1,:);
cnn_MSE_pos_ave(:, numel(x_control) + 1) = ...
    cnn_MSE_pos_ave(:,1);
cnn2_g2_MSE_pos_ave(numel(y_control) + 1, :) = ...
    cnn2_g2_MSE_pos_ave(1,:);
cnn2_g2_MSE_pos_ave(:, numel(x_control) + 1) = ...
    cnn2_g2_MSE_pos_ave(:,1);
cnn_g2_MSE_pos_ave(numel(y_control) + 1, :) = ...
    cnn_g2_MSE_pos_ave(1,:);
cnn_g2_MSE_pos_ave(:, numel(x_control) + 1) = ...
    cnn_g2_MSE_pos_ave(:,1);

% ***** Correlation_pos *****
figure(1)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
% set alpha values
alphaVal = 0.6;
alphaVal_colorbar = 0.8;

subplot(1,3,1)
ss = surf(x_control_surf, y_control_surf, real(cnn_correlation_pos_ave),...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title(' Average Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('CNN', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_corr);
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

subplot(1,3,2)
ss = surf(x_control_surf, y_control_surf, real(cnn_g2_correlation_pos_ave),...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Average Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('Group. CNN', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_corr);
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

subplot(1,3,3)
ss = surf(x_control_surf, y_control_surf, real(cnn2_g2_correlation_pos_ave),...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Average Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('Group. CNN 2x', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_corr);
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

% saving
saveas(figure(1), ...
    'Correlation_pos.png');

% ***** MSE_pos *****
figure(2)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.

subplot(1,3,1)
ss = surf(x_control_surf, y_control_surf, cnn_MSE_pos_ave,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title(' Average MSE', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('CNN', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_MSE);
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


subplot(1,3,2)
ss = surf(x_control_surf, y_control_surf, cnn_g2_MSE_pos_ave,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Average MSE', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('Group. CNN', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_MSE);
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


subplot(1,3,3)
ss = surf(x_control_surf, y_control_surf, cnn2_g2_MSE_pos_ave,...
    'FaceAlpha', alphaVal);
ss.EdgeColor = 'none';
title('Average MSE', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('Group. CNN 2x', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('x [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('y [m]', 'Interpreter', 'latex', 'FontSize', 15);
axis equal;
caxis(min_max_MSE);
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


% saving
saveas(figure(2), ...
    'MSE_pos.png');

% ***** Correlation_f *****
figure(3)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.

subplot(1,2,1)
% cnn1
plot(f, real(cnn_correlation_f_ave),'m*-', ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90],'MarkerSize', 8, 'linewidth', 1.5);
title('Average and Median Correlation', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('{Re\{.\}}', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, cnn_correlation_f_med, 'm^', 'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn_correlation_f_trm, 'm+', 'MarkerSize', 10, 'linewidth', 2)

% cnn group 2
plot(f, real(cnn_g2_correlation_f_ave),'*-', 'Color', [0.9290 0.4940 0.0950], ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
plot(f, cnn_g2_correlation_f_med, '^', 'Color', [0.9290 0.4940 0.0950], ...
    'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn2_correlation_f_trm, '+', 'Color', [0.9290 0.4940 0.0950], ...
%     'MarkerSize', 10, 'linewidth', 2)

% cnn2 group 2
plot(f, real(cnn2_g2_correlation_f_ave),'k*-', 'MarkerEdgeColor', ...
    [0.20, 0.75, 0.90], 'MarkerSize', 8, 'linewidth', 1.5);
plot(f, cnn2_g2_correlation_f_med, 'k^', 'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn_g2_correlation_f_trm, 'k+', 'MarkerSize', 10, 'linewidth', 2)
hold off
legend(['CNN Average. Global Average = '...
    num2str(real(cnn_correlation_f_glob_ave))], 'CNN Median', ...
    ['Group. CNN Average. Global Average = '...
    num2str(real(cnn_g2_correlation_f_glob_ave))], ...
    'Group. CNN Median', ['Group. CNN 2x Average. Global Average = ' ...
    num2str(real(cnn2_g2_correlation_f_glob_ave))], ...
    'Group. CNN 2x Median','location', 'SouthEast', 'Interpreter', ...
    'latex', 'FontSize', 12)

% ***** MSE_f *****
subplot(1,2,2)
% cnn1
plot(f, cnn_MSE_f_ave,'m*-', ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90],'MarkerSize', 8, 'linewidth', 1.5);
title('Average and Median MSE', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16, 'Interpreter', 'latex')
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('MSE [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, cnn_MSE_f_med, 'm^', 'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn_correlation_f_trm, 'm+', 'MarkerSize', 10, 'linewidth', 2)

% cnn group 2
plot(f, cnn_g2_MSE_f_ave,'*-', 'Color', [0.9290 0.4940 0.0950], ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
plot(f, cnn_g2_MSE_f_med, '^', 'Color', [0.9290 0.4940 0.0950], ...
    'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn2_correlation_f_trm, '+', 'Color', [0.9290 0.4940 0.0950], ...
%     'MarkerSize', 10, 'linewidth', 2)

% cnn2 group 2
plot(f, cnn2_g2_MSE_f_ave,'k*-', 'MarkerEdgeColor', [0.20, 0.75, 0.90], ...
    'MarkerSize', 8, 'linewidth', 1.5);
plot(f, cnn2_g2_MSE_f_med, 'k^', 'MarkerSize', 10, 'linewidth', 2)
% plot(f, cnn_g2_correlation_f_trm, 'k+', 'MarkerSize', 10, 'linewidth', 2)
hold off
legend(['CNN Average. Global Average = '...
    num2str(real(cnn_MSE_f_glob_ave)) ' dB'], 'CNN Median', ...
    ['Group. CNN Average. Global Average = '...
    num2str(real(cnn_g2_MSE_f_glob_ave)) ' dB'], ...
    'Group. CNN Median', ['Group. CNN 2x Average. Global Average = ' ...
    num2str(real(cnn2_g2_MSE_f_glob_ave)) ' dB'], ...
    'Group. CNN 2x Median','location', 'SouthEast', 'Interpreter', ...
    'latex', 'FontSize', 12)

% saving
saveas(figure(3), 'Correlation_f.png');
