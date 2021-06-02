% ***** ANC EVALUATION *****

% ***** Generate source and control region grids *****
% target field source- x-position variable
x_main_sources = -80:10:-20;
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

% ***** LOAD DATA *****
% ***** all-norm net *****
% init 
AE_allnorm = zeros(f_bins, numel(x_main_sources));
AE_allnorm_reg = zeros(f_bins, numel(x_main_sources));
ANC_power_allnorm = zeros(f_bins, numel(x_main_sources));
ANC_power_allnorm_reg = zeros(f_bins, numel(x_main_sources));
IL_allnorm = zeros(f_bins, numel(x_main_sources));
IL_allnorm_reg = zeros(f_bins, numel(x_main_sources));
main_power_allnorm = zeros(f_bins, numel(x_main_sources));
% load data
for i = 1:numel(x_main_sources)
    % change folder path. Files are stored in the "Distance" folder
    AE_allnorm(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/AE ' num2str(x_main_sources(i)) '.mat']);
    AE_allnorm_reg(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/AE_reg ' num2str(x_main_sources(i)) '.mat']);
    ANC_power_allnorm(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/ANC power ' num2str(x_main_sources(i)) '.mat']);
    ANC_power_allnorm_reg(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/ANC_reg power ' num2str(x_main_sources(i)) '.mat']);
    IL_allnorm(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/IL ' num2str(x_main_sources(i)) '.mat']);
    IL_allnorm_reg(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/IL_reg ' num2str(x_main_sources(i)) '.mat']);
    main_power_allnorm(:,i) = importdata(['Distance/' ...
        num2str(x_main_sources(i)) '/main power ' num2str(x_main_sources(i)) '.mat']);
end

% ***** Compute Averages *****
% need to convert back from log and the recompute log
AE_allnorm_ave_f = 10*log10(sum(10.^(AE_allnorm/10), 2)/numel(x_main_sources));
AE_allnorm_ave_d = 10*log10(sum(10.^(AE_allnorm/10), 1)/f_bins);
AE_allnorm_ave = 10*log10(sum(10.^(AE_allnorm_ave_f/10))/f_bins);

AE_allnorm_reg_ave_f = 10*log10(sum(10.^(AE_allnorm_reg/10), 2)/numel(x_main_sources));
AE_allnorm_reg_ave_d = 10*log10(sum(10.^(AE_allnorm_reg/10), 1)/f_bins);
AE_allnorm_reg_ave = 10*log10(sum(10.^(AE_allnorm_reg_ave_f/10))/f_bins);

ANC_power_allnorm_ave_f = 10*log10(sum(10.^(ANC_power_allnorm/10), 2)/...
    numel(x_main_sources));
ANC_power_allnorm_ave_d = 10*log10(sum(10.^(ANC_power_allnorm/10), 1)/f_bins);
ANC_power_allnorm_ave = 10*log10(sum(10.^(ANC_power_allnorm_ave_f/10))/f_bins);

ANC_power_allnorm_reg_ave_f = 10*log10(sum(10.^(ANC_power_allnorm_reg/10), 2)/...
    numel(x_main_sources));
ANC_power_allnorm_reg_ave_d = 10*log10(sum(10.^(ANC_power_allnorm_reg/10),...
    1)/f_bins);
ANC_power_allnorm_reg_ave = 10*log10(sum(10.^(ANC_power_allnorm_reg_ave_f/10))/f_bins);

IL_allnorm_ave_f = 10*log10(sum(10.^(IL_allnorm/10), 2)/numel(x_main_sources));
IL_allnorm_ave_d = 10*log10(sum(10.^(IL_allnorm/10), 1)/f_bins);
IL_allnorm_ave = 10*log10(sum(10.^(IL_allnorm_ave_f/10))/f_bins);

IL_allnorm_reg_ave_f = 10*log10(sum(10.^(IL_allnorm_reg/10), 2)/numel(x_main_sources));
IL_allnorm_reg_ave_d = 10*log10(sum(10.^(IL_allnorm_reg/10), 1)/f_bins);
IL_allnorm_reg_ave = 10*log10(sum(10.^(IL_allnorm_reg_ave_f/10))/f_bins);

main_power_allnorm_ave_f = 10*log10(sum(10.^(main_power_allnorm/10), 2)/...
    numel(x_main_sources));
main_power_allnorm_ave_d = 10*log10(sum(10.^(main_power_allnorm/10), 1)/f_bins);
main_power_allnorm_ave = 10*log10(sum(10.^(main_power_allnorm_ave_f/10))/f_bins);

% ***** PLOT *****
figure(1)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.
% marker styles
linM = {'o','+','*', '.', 'x', '_', '|'};

% IL
subplot(1,3,1)
plot(f, IL_allnorm_ave_f,'m-', 'linewidth', 6);
title('IL [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, IL_allnorm_reg_ave_f,'k-', 'linewidth', 6);
p1 = plot(f, IL_allnorm(:,1),'mo-', 'linewidth', 1.5);
% transparency 
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, IL_allnorm(:,i),'m', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
p1 = plot(f, IL_allnorm_reg(:,1),'ko-', 'linewidth', 1.5);
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, IL_allnorm_reg(:,i),'k', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
hold off
legend(['CNN Array. Average IL = ' num2str(IL_allnorm_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average IL = ' num2str(IL_allnorm_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex'); 

% AE
subplot(1,3,2)
plot(f, AE_allnorm_ave_f,'m-', 'linewidth', 6);
title('AE [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, AE_allnorm_reg_ave_f,'k-', 'linewidth', 6);
p1 = plot(f, AE_allnorm(:,1),'mo-', 'linewidth', 1.5);
% transparency 
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, AE_allnorm(:,i),'m', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
p1 = plot(f, AE_allnorm_reg(:,1),'ko-', 'linewidth', 1.5);
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, AE_allnorm_reg(:,i),'k', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
hold off
legend(['CNN Array. Average AE = ' num2str(AE_allnorm_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average AE = ' num2str(AE_allnorm_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex');

% Power
subplot(1,3,3)
plot(f, main_power_allnorm_ave_f,'-', 'Color', [0.20, 0.75, 0.90], ...
    'linewidth', 6);
title('Sound Field Power', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Power Level [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(f, ANC_power_allnorm_ave_f,'m-', 'linewidth', 6);
plot(f, ANC_power_allnorm_reg_ave_f,'k-', 'linewidth', 6);
p1 = plot(f, main_power_allnorm(:,1),'o-', 'linewidth', 1.5, ...
    'Color', [0.20, 0.75, 0.90]);
% transparency 
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, main_power_allnorm(:,i), 'linewidth', 1.5, 'Marker',linM{i}, ...
       'Color', [0.20, 0.75, 0.90]);
   p1.Color(4) = 0.25;
end
p1 = plot(f, ANC_power_allnorm(:,1),'mo-', 'linewidth', 1.5);
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, ANC_power_allnorm(:,i),'m', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
p1 = plot(f, ANC_power_allnorm_reg(:,1),'ko-', 'linewidth', 1.5);
p1.Color(4) = 0.25;
for i = 2:numel(x_main_sources)
   p1 = plot(f, ANC_power_allnorm_reg(:,i),'k', 'linewidth', 1.5, 'Marker',linM{i});
   p1.Color(4) = 0.25;
end
hold off
legend(['Main Array. Average Power = ' num2str(main_power_allnorm_ave) 'dB'], ...
    ['ANC CNN Array. Average Power = ' num2str(ANC_power_allnorm_ave) 'dB'], ...
    ['ANC Uni. Array. Average Power = ' num2str(ANC_power_allnorm_reg_ave) 'dB'], ...
    'location', 'NorthEast', 'FontSize', 11, 'Interpreter', 'latex');
% saving
saveas(figure(1), 'freq.png');

figure(2)
% figure bigger
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');           % set bg color to white.

% IL
subplot(1,3,1)
plot(x_main_sources, IL_allnorm_ave_d,'m*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
title('IL [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Main Array Location [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(x_main_sources, IL_allnorm_reg_ave_d,'k*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
hold off
legend(['CNN Array. Average AE = ' num2str(AE_allnorm_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average AE = ' num2str(AE_allnorm_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex');

% AE
subplot(1,3,2)
plot(x_main_sources, AE_allnorm_ave_d,'m*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
title('AE [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Main Array Location [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(x_main_sources, AE_allnorm_reg_ave_d,'k*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
hold off
legend(['CNN Array. Average AE = ' num2str(AE_allnorm_ave) 'dB'], ...
    ['Uniformly Spaced Array. Average AE = ' num2str(AE_allnorm_reg_ave) 'dB'], ...
    'location', 'SouthEast', 'FontSize', 11, 'Interpreter', 'latex');

% Power
subplot(1,3,3)
plot(x_main_sources, main_power_allnorm_ave_d,'*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8, ...
    'Color', [0.20, 0.75, 0.90]);
title('Sound Field Power [dB]', 'FontSize', 20, 'FontWeight', 'normal', ...
    'Interpreter', 'latex'); 
subtitle('', 'FontSize', 16)
xlabel('Main Array Location [m]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Mag [dB]', 'Interpreter', 'latex', 'FontSize', 15);
grid on;
hold on
plot(x_main_sources, ANC_power_allnorm_ave_d,'m*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
plot(x_main_sources, ANC_power_allnorm_reg_ave_d,'k*-', 'linewidth', 1.5, ...
    'MarkerEdgeColor', [0.20, 0.75, 0.90], 'MarkerSize', 8);
hold off
legend(['Main Array. Average Power = ' num2str(main_power_allnorm_ave) 'dB'], ...
    ['ANC CNN. Array. Average Power = ' num2str(ANC_power_allnorm_ave) 'dB'], ...
    ['ANC Uni. Array. Average Power = ' num2str(ANC_power_allnorm_reg_ave) 'dB'], ...
    'location', 'NorthEast', 'FontSize', 11, 'Interpreter', 'latex');
% saving
saveas(figure(2), 'dist.png');
