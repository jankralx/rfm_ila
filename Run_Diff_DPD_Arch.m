global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_ph c_pa_wp;

% PA coefficients selection
switch c_pa_sel
    case 1
        PA_Model_func = @(tx_sig, coefs, P, M) (MP_Output(tx_sig, P, M, coefs));
        % memoryless NXP_BLF8888A_75W
        c_pa_coefs = [0.225940289057156 - 1.30309938687855i;-0.753241844225673 + 1.50197079361271i;7.56239923012568 - 13.7974926416052i;-27.0418734124287 + 41.1956117815024i;44.6670818656791 - 53.4230010949243i;-35.0820341142530 + 32.0411865968688i;10.5932691363727 - 7.19328579093337i];
        c_pa_P = 7;
        c_pa_M = 0;
        c_pa_wp = 1;
        c_pa_gain = 1;
        c_pa_ph = 0.82;

    case 2
        PA_Model_func = @(tx_sig, coefs, P, M) (MP_Output(tx_sig, P, M, coefs));
        % memory NXP_BLF8888A_75W
        c_pa_coefs=[0.220916639043604 - 1.22044839893226i;-0.804479055683704 + 1.60679189505473i;7.32058974039844 - 14.7631487182526i;-24.7480079102736 + 41.5622939434676i;39.3595861467446 - 50.8984460940822i;-30.1038317346279 + 28.6281890757808i;8.91227214933715 - 5.90920863525877i;0.00563026755230095 - 0.0796768263048817i;0.0401168058092848 - 0.279616691127796i;0.447121832088464 + 1.99835120882601i;-3.35762879816203 - 2.76983932779456i;7.56741517081975 + 0.179412360318367i;-7.10434722256177 + 1.94725336144942i;2.41728650583326 - 0.978825330397617i];
        c_pa_P=7;
        c_pa_M=1;
        c_pa_wp = 1;
        c_pa_gain = 1;
        c_pa_ph = 0.8;

    case 3
        PA_Model_func = @(tx_sig, coefs, P, M) (MP_Output(tx_sig, P, M, coefs));
        % memoryless USRP RACOM PA
        c_pa_coefs=[0.845654964236630 + 1.49110374573613i;-2.76262037337481 + 1.88247060093166i;17.0635003318257 - 11.9190158113625i;-51.8235490530484 + 29.6408107759511i;81.3986654536564 - 40.1728308743692i;-64.9809973172151 + 27.2732197670390i;20.6182249621985 - 7.32252535080311i];
        c_pa_P=7;
        c_pa_M=0;
        c_pa_wp = 1;
        c_pa_gain = 1;
        c_pa_ph = 0.8;
        
    case 4
        PA_Model_func = @(tx_sig, coefs, P, M) (MP_Output(tx_sig, P, M, coefs));
        % memory USRP RACOM PA
        c_pa_coefs=[-1.28358468551180 + 0.0843179954454246i;3.28156363811602 + 0.538861754784714i;-51.9497317246123 - 10.8758256774022i;197.296500099173 + 63.8258872940854i;-325.030807619629 - 138.760958030307i;251.683328568101 + 129.840197319964i;-74.9973083479402 - 44.4054962278095i;-0.211979740065843 - 0.254061316215991i;2.81265481573720 + 1.05294498907392i;-16.9291176517167 - 15.8032198796011i;49.1238161078635 + 69.4126537921216i;-72.2393131692172 - 130.482761132544i;52.1519565108815 + 112.589070861896i;-14.7130368058352 - 36.6762110473874i];
        c_pa_P=7;
        c_pa_M=1;
        c_pa_wp = 1;
        c_pa_gain = 1;
        c_pa_ph = 3;
        
    otherwise
        PA_Model_func = @(tx_sig, coefs, P, M) (PA_Model(tx_sig.', coefs, P, M).');
        c_pa_coefs = [0.509360227527429 + 1.09524853022532i; -0.0992873613403637 - 0.170261849774516i; -0.0347754375003473 - 0.0247212149015436i; -0.00353320874772281 - 0.00211119148781448i; 0.00260430842062743 - 0.00429101487393531i; 0.00320810224865987 - 0.000580829859014498i; -0.000816817963483357 + 0.000357784194921971i];
        c_pa_P = 7;
        c_pa_M = 0;
        c_pa_wp = 0.8;
        c_pa_gain = 1.2;
        c_pa_ph = -1.13;
end

% estimate PA gain
if false
    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, c_mod_params);
    else
        tx_signal = FBMC_modulator(c_mod_params.num_subchannels, ...
            c_mod_params.oversample_fact, c_mod_params.num_frames);
    end
    tx_signal = tx_signal * c_pa_wp;
    
    % model PA for non-DPD output
    tx_signal_pa = PA_Model_func(tx_signal, c_pa_coefs, c_pa_P, c_pa_M);
    gain = lscov(tx_signal, tx_signal_pa);
    c_pa_gain = abs(gain);

    %Tx_pwr = sum(abs(tx_signal))
    %Tx_PA_pwr = sum(abs(tx_signal_pa))
    
    % now amplify the tx_signal by absolute value of gain
    %tx_signal = tx_signal * c_pa_gain;
    %c_pa_gain = 1;
end

% create tx signal -  we need dimensions of frequency axis
if strcmp(c_mod_params.type, 'QAM')
    binin = rand(c_mod_params.num_bits, 1) > 0.5;
    tx_signal = QAM_Modulator(binin, c_mod_params);
else
    tx_signal = FBMC_modulator(c_mod_params.num_subchannels, ...
        c_mod_params.oversample_fact, c_mod_params.num_frames);
end

freq_axis = linspace(-0.5, 0.5, length(tx_signal)+1);
freq_axis = freq_axis(2:end);


if ~c_use_indirect_learning_coefs
    % first coefficients constant
    dpd_num_coefs = size(dpd_par.Matrix_Func((1:10).', dpd_par.P, dpd_par.M), 2);
    coefs = zeros(dpd_num_coefs,1);
    coefs(1) = 1;
else
    % first coefficients estimated by indirect learning method for first time
    
    % firstly estimate coefficients using indirect learning (by the article
    % this is done only once in laboratory or derived from the AM/AM, AM/PM
    % characteristics).

    % create tx signal for transmission
    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, c_mod_params);
    else
        tx_signal = FBMC_modulator(c_mod_params.num_subchannels, ...
            c_mod_params.oversample_fact, c_mod_params.num_frames);
    end
    tx_signal = tx_signal/max(abs(tx_signal));      % ensure that the maximum value in tx signal is 1
    
    % model PA for non-DPD output
    tx_signal_pa = PA_Model_func(tx_signal*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    
    % calculate DPD coefficients
    U = DDR2_Matrix(tx_signal_pa, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal);
    
    % add some noise to coefficients to simulate time change
    coefs = coefs + abs(coefs).*(rand(size(coefs))-0.5)*0.1;
end

arch_names = {
    'ILA';
    'DLA';
    'R-DLA';
    'FM-ILA';
    'R-FM-ILA';
    'No DPD';
    };
arch_funcs = {
    @(tx_signal, coefs, dpd_par)(ila_dpd(tx_signal, coefs, dpd_par));
    @(tx_signal, coefs, dpd_par)(dla_dpd(tx_signal, coefs, dpd_par));
    @(tx_signal, coefs, dpd_par)(rdla_dpd(tx_signal, coefs, dpd_par));
    @(tx_signal, coefs, dpd_par)(fmila_dpd(tx_signal, coefs, dpd_par));
    @(tx_signal, coefs, dpd_par)(oqfmila_dpd(tx_signal, coefs, dpd_par));
    @(tx_signal, coefs, dpd_par)(nodpd(tx_signal, coefs));
    };

def.name = '';
def.func = @(tx_signal, coefs, dpd_par)(nodpd(tx_signal));
def.nmse = zeros(c_num_iter,1);
def.coefs = zeros(c_num_iter+1,length(coefs));
def.coefs(1,:) = coefs;
def.fftPA_sig_max = 0;
def.fftPA_sig_avg = 0;
def.mchan_pwr = zeros(c_num_iter,1);
def.acpr = zeros(c_num_iter,2,2);

% set starting mu
dpd_par.mu = c_mu;

% create parameter array from the default values
arch = repmat(def,length(arch_names),1);
for i = 1:size(arch,1)
    arch(i).name = arch_names{i};
    arch(i).func = arch_funcs{i};
    
    if strcmp(arch_names{i}, 'DLA') || strcmp(arch_names{i}, 'R-DLA')
        arch(i).coefs(1,1) = 0.2;
    end        
end




for iter = 1:c_num_iter
    % create tx signal for transmission
    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, c_mod_params);
    else
        tx_signal = FBMC_modulator(c_mod_params.num_subchannels, ...
            c_mod_params.oversample_fact, c_mod_params.num_frames);
    end
    tx_signal = tx_signal/max(abs(tx_signal));      % ensure that the maximum value in tx signal is 1
    
    for i = 1:size(arch,1)
        coefs = arch(i).coefs(iter,:).';
        [sig_pa, err, coefs] = arch(i).func(tx_signal, coefs, dpd_par);

        % calculate ACPRs and power in the main channel
        [acpr_val, mchan_pwr] = ACPR(sig_pa, ...
            c_mod_params.oversample_fact, c_acrp_channels, false);

        if i == length(arch_names)   % backoff only for PA without DPD
            backoff = 10^((mchan_pwr-arch(4).mchan_pwr(iter))/20);
            [sig_pa, err, coefs] = arch(i).func(tx_signal/backoff, coefs, dpd_par);
        
            % calculate ACPRs and power in the main channel
            [acpr_val, mchan_pwr] = ACPR(sig_pa, ...
                c_mod_params.oversample_fact, c_acrp_channels, false);
        end
        
        % save parameters
        arch(i).acpr(iter,:,:) = acpr_val;
        arch(i).mchan_pwr(iter) = mchan_pwr;

        arch(i).nmse(iter) = err;
        arch(i).coefs(iter+1,:) = coefs.';

        if c_is_plot || iter > c_avg_fft_after_iter
            Sig_pa = fft(sig_pa);       % calculate spectrum if it is needed
        end
    
        if iter > c_avg_fft_after_iter
            arch(i).fftPA_sig_avg = arch(i).fftPA_sig_avg + abs(Sig_pa);
            arch(i).fftPA_sig_max = max(arch(i).fftPA_sig_max.', abs(Sig_pa).').';
        end
        
        if c_is_plot
            figure(1);
            if i > 1
                hold on;
            end
            plot(freq_axis, db(fftshift(Sig_pa)), 'DisplayName',arch(i).name);
            hold off;

            % NMSE evolution
            if c_num_iter > 1
                figure(2);
                if i > 1
                    hold on;
                end
                plot(arch(i).nmse(1:iter), '-x', 'DisplayName',arch(i).name);
                hold off;
            else
                figure(1);
                text(0, 20-i*5, sprintf(' %s NMSE = %.1f', arch(i).name, arch(i).nmse(1)));
            end;
            
            % ACPR evolution
            if c_num_iter > 1
                figure(3);
                if i == 1
                    clf;
                end
                subplot(2,1,1);
                hold on;
                plot(Avg_dB(arch(i).acpr(1:iter,1,:),3,10), '-x', 'DisplayName',arch(i).name);
                hold off;

                subplot(2,1,2);
                hold on;
                plot(Avg_dB(arch(i).acpr(1:iter,2,:),3,10), '-x', 'DisplayName',arch(i).name);
                hold off;
            end            
        end
        
        % decrease the step size
        if iter > c_mu_change_after
            dpd_par.mu = dpd_par.mu / c_mu_gamma;
        end
    end
            
   
    % set plot axis and titles
    figure(1);
    % move no DPD amplifier into bottom
    hLines = get(gca,'Children');
    temp = hLines(1);
    hLines(1) = hLines(end);
    hLines(end) = temp;
    set(gca,'Children',hLines)
    legend('show');
    title('Magnitude spectra');
    xlabel('Frequency');
    ylabel('Magnitude (dB)');
    axis([-0.5 0.5, -100 60]);
    
    if c_num_iter > 1
        figure(2);
        legend('show');
        title('Normalised mean square error');
        xlabel('Iteration cycle');
        ylabel('NMSE (dB)');
        axis([0 c_num_iter, -60 0]);
    end;
    
    if c_num_iter > 1
        figure(3);
        subplot(2,1,1);
        legend('show');
        title('ACPR in 1-st adjacent channel');
        xlabel('Iteration cycle');
        ylabel('ACPR (dB)');
        axis([0 c_num_iter, -70 0]);
        subplot(2,1,2);
        legend('show');
        title('ACPR in 2-nd adjacent channel');
        xlabel('Iteration cycle');
        ylabel('ACPR (dB)');
        axis([0 c_num_iter, -70 0]);
    end;
    
    drawnow;
end

% average saved spectra
for i = 1:size(arch,1)
    arch(i).fftPA_sig_avg = arch(i).fftPA_sig_avg / (c_num_iter - c_avg_fft_after_iter);
end

% plot result spectra
figure(1);
clf;
for i = 1:size(arch,1)
    hold on;%     subplot(1,2,2);
    plot(freq_axis, db(fftshift(arch(i).fftPA_sig_avg)), 'DisplayName',arch(i).name);
    hold off;
end

legend('show');
title('Magnitude spectra');
xlabel('Frequency');
ylabel('Magnitude (dB)');
axis([-inf inf, -30 60]);

%% save plots
if (1)
    am_fig = figure(5);
    am_fig.Renderer='Painters';
    set(gcf, 'Position', [0 0 600 400]);
    title('\textbf{AM/AM Characteristics}');
    xlabel('Relative Input Level (-)');
    ylabel('Relative Output Level (-)');
    axis([0 1 0 1]);
    grid on;
    
%    annotation('textarrow',[0.40 0.43],[0.7 0.66],'String','PA model');
%     annotation('textarrow',[0.28 0.31],[0.53 0.48],'String','Real PA');
%     annotation('textarrow',[0.535 0.55],[0.63 0.57],'String','Linearised');
%     annotation('textarrow',[0.66 0.675],[0.53 0.47],'String','DPD output');
    
    ApplyFigureSettings(am_fig);
    saveas(gcf, 'figures/am_am.pdf');
end


%%


% iterative learning control
function [tx_signal_pa, nmserr] = oqilc_dpd(tx_signal, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_wp c_pa_gain;

    c_iter_lim = 300;
    % iterative learning control loop
    ilc_nmse = [];
    tx_dpd = tx_signal;
    tx_signal_pa = PA_Model_func(tx_dpd/10, c_pa_coefs, c_pa_P, c_pa_M);
    gss = sum(abs(tx_signal_pa))/sum(abs(tx_dpd)/10);
    gamma = 2/gss
    
    h = waitbar(0,'Please wait...');
    for i = 1:c_iter_lim
        tx_signal_pa = PA_Model_func(tx_dpd, c_pa_coefs, c_pa_P, c_pa_M);
        tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;

        % linear
        tx_dpd = tx_dpd + gamma*(tx_signal-tx_signal_pa);
        
        % gained-based
%         g = tx_dpd./tx_signal_pa;
%         g(isnan(g)) = 0;
%         tx_dpd = tx_dpd + (tx_signal-tx_signal_pa).*g/5;
                
%         tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
        ilc_nmse(i) = nmse(tx_signal, tx_signal_pa);
        if ilc_nmse(i) < -40 && ilc_nmse(i-1) < ilc_nmse(i)
            break;
        end
        waitbar(i / c_iter_lim);
    end
    
    % plot evalulation of the ILC NMSE
    figure(4);
    plot(ilc_nmse);
    
    close(h);
    
    nmserr = ilc_nmse(end);
end

% iterative learning control
function [tx_signal_pa, nmserr] = ilc_dpd(tx_signal, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_wp c_pa_gain c_mod_params;

    c_iter_lim = 300;

    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_sig_t = QAM_Modulator(binin, c_mod_params);
    else
        tx_sig_t = FBMC_modulator(c_mod_params.num_subchannels, ...
            c_mod_params.oversample_fact, c_mod_params.num_frames);
    end
    
    tx_sig_t = tx_sig_t / c_pa_wp;

    ilc_nmse = [];
    tx_dpd = tx_sig_t;
%     tx_signal_pa = PA_Model_func((tx_dpd*c_pa_wp/10).', c_pa_coefs, c_pa_P, c_pa_M).';
%     gss = sum(abs(tx_signal_pa))/sum(abs(tx_dpd)/10);
%     gamma = 2/gss;
    
    % iterative learning control loop    
    h = waitbar(0,'Please wait...');
    for i = 1:c_iter_lim
        tx_signal_pa = PA_Model_func((tx_dpd).', c_pa_coefs, c_pa_P, c_pa_M).';
%         tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;

        % linear
%         tx_dpd = tx_dpd + gamma*(tx_signal-tx_signal_pa);
        
        % gained-based
        g = tx_dpd./tx_signal_pa;
        g(isnan(g)) = 0;
        tx_dpd = tx_dpd + (tx_sig_t-tx_signal_pa/c_pa_gain).*g/5;
                
%         tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
        ilc_nmse(i) = nmse(tx_sig_t, tx_signal_pa);
        if ilc_nmse(i) < -40 && ilc_nmse(i-1) < ilc_nmse(i)
            break;
        end
        waitbar(i / c_iter_lim);
    end
    
    % plot evalulation of the ILC NMSE
    figure(4);
    subplot(2,1,1);
    plot(ilc_nmse);
    subplot(2,1,2);
    plot(db(abs(fftshift(fft(tx_signal_pa)))));
    
    % from the trained PA input, calculate the model of the DPD
    U = DDR2_Matrix(tx_sig_t, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_dpd);
    
    % with the DPD model predistort the input tx_signal and evalueate it
    % on the PA
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model_func((dpd_out*c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
    
    close(h);
end

% forward PA model extraction and subsequent DPD calculation
function [tx_signal_pa, nmserr, coefs] = oqfmila_dpd(tx_signal, coefs, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_wp c_pa_gain c_params;
    
    % predistort transmitted signal with so-far known coefficients
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    dpd_out = U * coefs;
    
    % transmit the signal through the PA model
    tx_signal_pa = PA_Model_func(dpd_out*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    tx_signal_pan = awgn(tx_signal_pa, c_params.fb.SNR, db(0.5,'power')); % add noise to feedback signal
    
    % calculate forward PA model coefficients
    U = dpd_par.Matrix_Func(dpd_out, dpd_par.P, dpd_par.M);
    M = [real(U) -imag(U)];
    b = lscov(M, real(tx_signal_pan));
    pa_coefs = b(1:(size(b,1)/2)) + 1j*b((size(b,1)/2)+1:end);
    
    % calculate forward PA model output
    U = dpd_par.Matrix_Func(dpd_out, dpd_par.P, dpd_par.M);
    tx_signal_pa_mod = U * pa_coefs;
  
    % exclude values which are greater than 1
    dpd_out_sel = dpd_out(abs(tx_signal_pa_mod) < 1);
    tx_signal_pa_mod = tx_signal_pa_mod(abs(tx_signal_pa_mod) < 1);


    % we have noise-free model output, calculate the required coefficients now
    U = dpd_par.Matrix_Func(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, dpd_out_sel);

    % plot AM/AM characteristisc
    figure(5);
    plot(abs(dpd_out), abs(tx_signal_pa), 'o',  'DisplayName','Real PA');
    hold on;
    plot(abs(dpd_out_sel), abs(tx_signal_pa_mod), '.',  'DisplayName','PA model');
    plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
    plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
    hold off;
    title('AM/AM characteristics');
    legend('show', 'Location','southeast');
       
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end


% forward PA model extraction and subsequent DPD calculation
function [tx_signal_pa, nmserr, coefs] = fmila_dpd(tx_signal, coefs, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_wp c_pa_gain c_pa_ph c_params;
    
    % predistort transmitted signal with so-far known coefficients
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    dpd_out = U * coefs;
    
    % transmit the signal through the PA model
    tx_signal_pa = PA_Model_func(dpd_out*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    tx_signal_pan = awgn(tx_signal_pa, c_params.fb.SNR, db(0.5,'power')); % add noise to feedback signal
    %tx_signal_pan(abs(tx_signal_pan)>1) = exp(1j*angle(tx_signal_pan(abs(tx_signal_pan)>1)));

    % calculate forward PA model coefficients
    U = dpd_par.Matrix_Func(dpd_out, dpd_par.P, dpd_par.M);
    pa_coefs  = lscov(U, tx_signal_pan);
        
    % calculate forward PA model output
    U = dpd_par.Matrix_Func(dpd_out, dpd_par.P, dpd_par.M);
    tx_signal_pa_mod = U * pa_coefs;
    
    % exclude values which are greater than 1
    dpd_out_sel = dpd_out(abs(tx_signal_pa_mod) < 1);
    tx_signal_pa_mod = tx_signal_pa_mod(abs(tx_signal_pa_mod) < 1);

    % we have noise-free model output, calculate the required coefficients now
    U = dpd_par.Matrix_Func(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, dpd_out_sel);
    
    % plot AM/AM characteristisc
%     figure(5);
% %     subplot(1,2,2);
%     plot(abs(dpd_out), abs(tx_signal_pa), '.',  'DisplayName','Real PA');
%     hold on;
%     plot(abs(dpd_out), abs(tx_signal_pa_mod), '.',  'DisplayName','DPD output');
%     plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
%     plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
%     hold off;
%     title('AM/AM characteristics');
%     legend('show', 'Location','southeast');
       
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end


function [tx_signal_pa, nmserr, coefs] = rdla_dpd(tx_signal, coefs, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_ph c_pa_wp c_params;

    % calculate DPD output
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    dpd_out = U * coefs;

    % calculate PA model output
    tx_signal_pa = PA_Model_func(dpd_out*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    tx_signal = tx_signal*exp(-1j*c_pa_ph); % phase needs to be approx. synchronised
    tx_signal_pan = awgn(tx_signal_pa, c_params.fb.SNR, db(0.5,'power')); % add noise to feedback signal

    % direct learning DPD with one quadrature
    % calculate model matrix
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    M = [real(U) -imag(U)];

    % use Gauss-Newton method to update coefficients
    b = [real(coefs); imag(coefs)];
    b = b - dpd_par.mu*((M'*M)\(M'*real(tx_signal_pan-tx_signal)));
    coefs = b(1:(size(b,1)/2)) + 1j*b((size(b,1)/2)+1:end);
        
    % plot AM/AM characteristisc
%     figure(4);
%     plot(abs(dpd_out), abs(tx_signal_pa), '.',  'DisplayName','PA');
%     hold on;
%     plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD');
%     plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
%     hold off;
%     title('AM/AM characteristics of the PA');
%     legend('show', 'Location','southeast');
    
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end


function [tx_signal_pa, nmserr, coefs] = dla_dpd(tx_signal, coefs, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_ph c_pa_wp c_params;

    % calculate DPD output
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    dpd_out = U * coefs;

    % calculate PA model output
    tx_signal_pa = PA_Model_func(dpd_out*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
%     g = lscov(tx_signal_pa, tx_signal)
%     an = angle(g)
%     amp = abs(g)
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    tx_signal_pa = exp(1j*(c_pa_ph)) * tx_signal_pa;   % phase needs to be approx. synchronised
    tx_signal_pan = awgn(tx_signal_pa, c_params.fb.SNR, db(0.5,'power')); % add noise to feedback signal
    
    % use Gauss-Newton method to update coefficients
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    coefs = coefs - dpd_par.mu*((U'*U)\(U'*(tx_signal_pan-tx_signal)));
        
    % plot AM/AM characteristisc
%     figure(3);
%     plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
%     hold on;
%     plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
%     hold off;
%     title('AM/AM characteristics');
%     legend('show', 'Location','southeast');
    
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end

function [tx_signal_pa, nmserr, coefs] = ila_dpd(tx_signal, coefs, dpd_par)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_wp c_params;

    % predistort transmitted signal with so-far known coefficients
    U = dpd_par.Matrix_Func(tx_signal, dpd_par.P, dpd_par.M);
    dpd_out = U * coefs;    
    
    % transmit the signal through the PA model
    tx_signal_pa = PA_Model_func(dpd_out*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    tx_signal_pan = awgn(tx_signal_pa, c_params.fb.SNR, db(0.5,'power')); % add noise to feedback signal
   
    % calculate new DPD coefficients
    U = dpd_par.Matrix_Func(tx_signal_pan, dpd_par.P, dpd_par.M);
    coefs = lscov(U, dpd_out);

    % plot AM/AM characteristisc
%     figure(4);
%     plot(abs(dpd_out), abs(tx_signal_pa), '.',  'DisplayName','PA');
%     hold on;    
%     plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
%     plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
%     hold off;
%     title('AM/AM characteristics');
%     legend('show', 'Location','southeast');   
       
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end

function [tx_signal_pa, nmserr, coefs] = nodpd(tx_signal, coefs)
    global PA_Model_func c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_wp;

    % model PA for non-DPD output
    tx_signal_pa = PA_Model_func(tx_signal*c_pa_wp, c_pa_coefs, c_pa_P, c_pa_M);
    tx_signal_pa = tx_signal_pa / c_pa_gain / c_pa_wp;
    
    % calculate NMSE
    tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pan);
end