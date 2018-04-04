addpath('Utils');
addpath('mutils');

%--------------------------------------------------------------------------
global c_mod_params c_params;

% set to true if plots are required
c_is_plot = true;

% Modulator selection
c_mod_params.type = 'FBMC';

% FBMC Modulator
c_mod_params.num_subchannels = 1024; % num of SUBCHANNELS
c_mod_params.oversample_fact = 6; % 12 % oversampling
c_mod_params.num_frames = 18; % 18

% QAM Modulator
c_mod_params.num_bits = 1024*6;
c_mod_params.M = 16;
c_mod_params.filt = 'sqrt';

% noise in feedback
c_params.fb.SNR = 30;

% Number of iterations
c_num_iter = 220;

c_avg_fft_after_iter = 20;

% PA model selection
c_pa_sel = 4;

c_acrp_channels = [1/2+0.1, 3/2+0.1,...
                   3/2+0.2, 5/2+0.2];

dpd_par.Matrix_Func = @(tx_sig, P, M)(MP_Matrix(tx_sig, P, M));
dpd_par.P = 7;
dpd_par.M = 3;

% direct learning convergence settings
c_mu = 0.8;                 % convergence speed
c_mu_change_after = 10;     % number of iterations after the mu is decreased
c_mu_gamma = 1.5;           % how many times the mu is decreased each iteration

c_use_indirect_learning_coefs = false;      % only for direct learning methods

Run_Diff_DPD_Arch();

% save results
arch = arch1;
save('results.mat', 'arch', 'freq_axis', 'c_avg_fft_after_iter');


%% change DPD parameters
% dpd_par.P = 3;
% dpd_par.M = 1;
% 
% Run_Diff_DPD_Arch();
% arch2 = arch;

%% sweep SNR to get the dependency of NMSE, ACPR on SNR
c_num_iter = 220;
c_avg_fft_after_iter = 20;

SNR_sweep = 0:2.5:65;
arch_SNR = {};

for k = 1:length(SNR_sweep)
    c_params.fb.SNR = SNR_sweep(k);
    fprintf('SNR = %i dB\n', c_params.fb.SNR);
    Run_Diff_DPD_Arch();
    arch_SNR{k} = arch;
end

% save('results_SNR.mat', 'arch_SNR', 'SNR_sweep', 'freq_axis', 'c_avg_fft_after_iter');


