close all;

% load results if they are not in the variables
if ~exist('arch', 'var')
    load('results.mat');
end

% define default colors  for plots
colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];



% plot average spectra
spect_fig = figure(1);
clf;
for i = 1:size(arch,1)
    hold on;
    plot(freq_axis, db(fftshift(arch(i).fftPA_sig_avg)), ...
        'DisplayName',arch(i).name, 'Color',colors(i,:));
    hold off;
end

% move no DPD amplifier into bottom
hLines = get(gca,'Children');
hLines = [hLines(2:end); hLines(1)];
set(gca,'Children',hLines);

% add text labels
annotation('textarrow',[0.337 0.365],[0.58 0.535],'String','No DPD');
annotation('textarrow',[0.32 0.30],[0.2 0.25],'String','ILA, DLA, R-DLA, FM-ILA, R-FM-ILA');

spect_fig.Renderer='Painters';
set(gcf, 'Position', [0 0 620 400]);
title('\textbf{Apmlitude Frequency Spectra}');
xlabel('Frequency ($F_S$)');
ylabel('Signal Magnitude (dB)');
axis([-0.5 0.5 -30 50]);
legend('show');
grid on;
ApplyFigureSettings(spect_fig);
saveas(gcf, 'figures/spectrum.pdf');

%%
% plot NMSE evolution
nmse_fig = figure(2);
clf;
marks = ['x','x','o','o','+','x'];
for i = 1:size(arch,1)
    hold on;
    plot(arch(i).nmse, sprintf('-%s',marks(i)), ...
        'DisplayName',arch(i).name, 'Color',colors(i,:));
    hold off;
end

% move no DPD amplifier into bottom
hLines = get(gca,'Children');
hLines = [hLines(2:end); hLines(1)];
set(gca,'Children',hLines);


% add text labels
annotation('textarrow',[0.28 0.25],[0.29 0.27],'String','ILA, FM-ILA, R-FM-ILA');
annotation('textarrow',[0.48 0.45],[0.6 0.57],'String','DLA, R-DLA');
annotation('textarrow',[0.48 0.45],[0.83 0.79],'String','No DPD');

spect_fig.Renderer='Painters';
set(gcf, 'Position', [0 0 600 400]);
title('\textbf{NMSE Evolution}');
xlabel('Iteration Cycle (-)');
ylabel('NMSE (dB)');
axis([0 15 -45 -15]);
legend('show', 'Location', 'East');
grid on;
ApplyFigureSettings(nmse_fig);
saveas(gcf, 'figures/nmse.pdf');


%% NMSE and ACPR in a table for latex
ilist = 1:size(arch,1);
ilist = [ilist(end), ilist(1:end-1)];
for i = ilist
    fprintf('%s & %.1f & %.1f & %.1f', arch(i).name, ...
        Avg_dB(arch(i).nmse(c_avg_fft_after_iter:end),1,10),...
        Avg_dB(Avg_dB(arch(i).acpr(c_avg_fft_after_iter:end,1,:),3,10),1,10),...
        Avg_dB(Avg_dB(arch(i).acpr(c_avg_fft_after_iter:end,2,:),3,10),1,10));
    fprintf(' \\\\\n\\hline \n');
end

ilist = 1:size(arch2,1);
ilist = [ilist(end), ilist(1:end-1)];
for i = ilist
    fprintf('%s & %.1f & %.1f & %.1f', arch(i).name, ...
        Avg_dB(arch2(i).nmse(c_avg_fft_after_iter:end),1,10),...
        Avg_dB(Avg_dB(arch2(i).acpr(c_avg_fft_after_iter:end,1,:),3,10),1,10),...
        Avg_dB(Avg_dB(arch2(i).acpr(c_avg_fft_after_iter:end,2,:),3,10),1,10));
    fprintf(' \\\\\n\\hline \n');
end


%% Noise influence

% calculate average NMSE and ACPR for all SNR settings

% nmse = zeros(size(arch_SNR,1),length(SNR_sweep));
% acpr_1 = zeros(size(arch_SNR,1),length(SNR_sweep));
% acpr_2 = zeros(size(arch_SNR,1),length(SNR_sweep));
% 
% for k = 1:length(SNR_sweep)
%     arch = arch_SNR{k};
%     for i = 1:size(arch,1)
%         nmse(i,k) = Avg_dB(arch(i).nmse(c_avg_fft_after_iter:end),1,10);
%         acpr_1(i,k) = Avg_dB(Avg_dB(arch(i).acpr(c_avg_fft_after_iter:end,1,:),3,10),1,10);
%         acpr_2(i,k) = Avg_dB(Avg_dB(arch(i).acpr(c_avg_fft_after_iter:end,2,:),3,10),1,10);
%     end
% end
% 
% figure(3);
% clf;
% 
% ilist = 1:size(arch,1);
% ilist = [ilist(end), ilist(1:end-1)];
% for i = ilist
%     hold on;
%     plot(SNR_sweep, nmse(i,:), 'DisplayName',arch(i).name, 'Color',colors(i,:));
%     hold off;
% end
% legend('show');






















