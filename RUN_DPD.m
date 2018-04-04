addpath('Utils');
addpath('mutils');

%--------------------------------------------------------------------------
global c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_wp c_c0 c_gamma c_comp_val c_is_prob_sel;

%FBMC Modulator
c_num_subchannels = 1024; % num of SUBCHANNELS
c_oversample_fact = 12; % 12 % oversampling
c_num_frames = 12; % 18

% QAM Modulator
c_mod_params.type = 'FBMC';
c_mod_params.num_bits = 1024*6;
c_mod_params.M = 16;
c_mod_params.filt = 'sqrt';
c_mod_params.oversample_fact = 12;

% Number of iterations
c_num_iter = 1;

% PA model selection
c_pa_sel = 1;

% switches on or off probability selection of points in AM/AM characteristics
c_is_prob_sel = false;

dpd_par.P = 5;
dpd_par.M = 0;
dpd_par.mu = 0.2;   % convergence speed

c_use_indirect_learning_coefs = false;      % only for direct learning methods

% input comparator levels
if strcmp(c_mod_params.type, 'QAM')
    c_comp_val = [0.33]; %-0.95:0.5:0.95;
else
    c_comp_val = [0.5]; %-0.95:0.5:0.95;
end;


% PA coefficients
switch c_pa_sel
    case 1
        % memoryless NXP_BLF8888A_75W
        c_pa_coefs = [0.225940289057156 - 1.30309938687855i;-0.753241844225673 + 1.50197079361271i;7.56239923012568 - 13.7974926416052i;-27.0418734124287 + 41.1956117815024i;44.6670818656791 - 53.4230010949243i;-35.0820341142530 + 32.0411865968688i;10.5932691363727 - 7.19328579093337i];
        c_pa_P = 7;
        c_pa_M = 0;
        c_pa_wp = 1/18;
        c_pa_gain = 1300;
        c_comp_val = 0.8;       % FBMC

    case 2
        % memory NXP_BLF8888A_75W
        c_pa_coefs=[0.220916639043604 - 1.22044839893226i;-0.804479055683704 + 1.60679189505473i;7.32058974039844 - 14.7631487182526i;-24.7480079102736 + 41.5622939434676i;39.3595861467446 - 50.8984460940822i;-30.1038317346279 + 28.6281890757808i;8.91227214933715 - 5.90920863525877i;0.00563026755230095 - 0.0796768263048817i;0.0401168058092848 - 0.279616691127796i;0.447121832088464 + 1.99835120882601i;-3.35762879816203 - 2.76983932779456i;7.56741517081975 + 0.179412360318367i;-7.10434722256177 + 1.94725336144942i;2.41728650583326 - 0.978825330397617i];
        c_pa_P=7;
        c_pa_M=1;
        c_pa_wp = 1/18;
        c_pa_gain = 1700;
        c_comp_val = 0.45;       % FBMC

    case 3
        % memoryless USRP RACOM PA
        c_pa_coefs=[0.845654964236630 + 1.49110374573613i;-2.76262037337481 + 1.88247060093166i;17.0635003318257 - 11.9190158113625i;-51.8235490530484 + 29.6408107759511i;81.3986654536564 - 40.1728308743692i;-64.9809973172151 + 27.2732197670390i;20.6182249621985 - 7.32252535080311i];
        c_pa_P=7;
        c_pa_M=0;
        c_pa_wp = 1/18;
        c_pa_gain = 2000;
        c_comp_val = 0.05;%[0.2 0.5 0.7 0.9]; %0.3;       % FBMC
        
    case 4
        % memory USRP RACOM PA
        c_pa_coefs=[-1.28358468551180 + 0.0843179954454246i;3.28156363811602 + 0.538861754784714i;-51.9497317246123 - 10.8758256774022i;197.296500099173 + 63.8258872940854i;-325.030807619629 - 138.760958030307i;251.683328568101 + 129.840197319964i;-74.9973083479402 - 44.4054962278095i;-0.211979740065843 - 0.254061316215991i;2.81265481573720 + 1.05294498907392i;-16.9291176517167 - 15.8032198796011i;49.1238161078635 + 69.4126537921216i;-72.2393131692172 - 130.482761132544i;52.1519565108815 + 112.589070861896i;-14.7130368058352 - 36.6762110473874i];
        c_pa_P=7;
        c_pa_M=1;
        c_pa_wp = 1/18;
        c_pa_gain = 9200;
        c_comp_val = 0.5;       % FBMC
    
    otherwise
        c_pa_coefs = [0.509360227527429 + 1.09524853022532i; -0.0992873613403637 - 0.170261849774516i; -0.0347754375003473 - 0.0247212149015436i; -0.00353320874772281 - 0.00211119148781448i; 0.00260430842062743 - 0.00429101487393531i; 0.00320810224865987 - 0.000580829859014498i; -0.000816817963483357 + 0.000357784194921971i];
        c_pa_P = 7;
        c_pa_M = 0;
        c_pa_wp = 0.9;
        c_pa_gain = 1;
        c_comp_val = 0.2;
end

% estimate PA gain
if false
    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, c_mod_params);
    else
        tx_signal = FBMC_modulator(c_num_subchannels, c_oversample_fact, c_num_frames);
    end
    tx_signal = tx_signal * c_pa_wp;
    
    % model PA for non-DPD output
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    gain = lscov(tx_signal, tx_signal_pa);
    c_pa_gain = abs(gain);

    %Tx_pwr = sum(abs(tx_signal))
    %Tx_PA_pwr = sum(abs(tx_signal_pa))
    
    % now amplify the tx_signal by absolute value of gain
    %tx_signal = tx_signal * c_pa_gain;
    %c_pa_gain = 1;
end


if ~c_use_indirect_learning_coefs
    % first coefficients constant
    dpd_num_coefs = floor((dpd_par.P-1)/2)*(4*dpd_par.M+1)+(dpd_par.M+1);
    coefs = zeros(dpd_num_coefs,1);
    coefs(1) = 0.5;
else
    % first coefficients estimated by indirect learning method for first time
    
    % firstly estimate coefficients using indirect learning (by the article
    % this is done only once in laboratory or derived from the AM/AM, AM/PM
    % characteristics).

    % create tx signal for transmission
    tx_signal = FBMC_modulator(c_num_subchannels, c_oversample_fact, c_num_frames);

    % calculate PA model output
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';

    coefs = calc_DPD_DDR2(tx_signal_pa, tx_signal, dpd_par.P, dpd_par.M);

    % add some noise to coefficients to simulate time change
    % coefs = coefs + abs(coefs).*rand(size(coefs))*0.2;
end

% -------------------------------------------------------------------------
% Direct learning approach 

% one bit DPD parameters
c_c0 = 0.07;
c_gamma = 0.9;

nodpd_nmserr = zeros(c_num_iter,1);

ilc_nmserr = zeros(c_num_iter,1);
uilc_nmserr = zeros(c_num_iter,1);
ruilc_nmserr = zeros(c_num_iter,1);
oq_ruilc_nmserr = zeros(c_num_iter,1);
oqc_ilc_nmserr = zeros(c_num_iter,1);
oq_ilc_nmserr = zeros(c_num_iter,1);

dl_nmserr = zeros(c_num_iter,1);
dl_coefs = coefs;

idl_nmserr = zeros(c_num_iter,1);

for iter = 1:c_num_iter
    % create tx signal for transmission
    if strcmp(c_mod_params.type, 'QAM')
        binin = rand(c_mod_params.num_bits, 1) > 0.5;
        tx_signal = QAM_Modulator(binin, c_mod_params);
    else
        tx_signal = FBMC_modulator(c_num_subchannels, c_oversample_fact, c_num_frames);
    end
    %tx_signal = tx_signal .* c_pa_wp;
    
    % model PA for non-DPD output
    tx_signal_pa = PA_Model((tx_signal*c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    err = nmse(tx_signal, tx_signal_pa);
    nodpd_nmserr(iter) = err;


    % different variants of DPD
    [oqc_ilc_sig_pa, err] = oqc_ilc_dpd(tx_signal, dpd_par);
    oqc_ilc_nmserr(iter) = err;

%     [oq_ruilc_sig_pa, err] = oq_ruilc_dpd(tx_signal, dpd_par);
%     oq_ruilc_nmserr(iter) = err;
% 
%     [oq_ilc_sig_pa, err] = oq_ilc_dpd(tx_signal, dpd_par);
%     oq_ilc_nmserr(iter) = err;

%     [ruilc_sig_pa, err] = ruilc_dpd(tx_signal, dpd_par);
%     ruilc_nmserr(iter) = err;
%     
%     [uilc_sig_pa, err] = uilc_dpd(tx_signal, dpd_par);
%     uilc_nmserr(iter) = err;
    
    [ilc_sig_pa, err] = ilc_dpd(tx_signal, dpd_par);
    ilc_nmserr(iter) = err;

%     [dl_sig_pa, err, dl_coefs] = dl_dpd(tx_signal, dl_coefs, dpd_par);
%     dl_nmserr(iter) = err;
    
%     [idl_sig_pa, err] = idl_dpd(tx_signal, dpd_par);
%     idl_nmserr(iter) = err;
    
    % plot results
    Tx_signal_pa = fft(tx_signal_pa);
    Oqc_ilc_sig_pa = fft(oqc_ilc_sig_pa);
%     Oq_ruilc_sig_pa = fft(oq_ruilc_sig_pa);
%     Oq_ilc_sig_pa = fft(oq_ilc_sig_pa);
%     Ruilc_sig_pa = fft(ruilc_sig_pa);
%     Uilc_sig_pa = fft(uilc_sig_pa);
    Ilc_sig_pa = fft(ilc_sig_pa);
%     Dl_sig_pa = fft(dl_sig_pa);
%     Idl_sig_pa = fft(idl_sig_pa);
    
    figure(1);
    % signals spectra
    if c_num_iter > 1
        subplot(2,1,1);
    else
        subplot(1,1,1);
    end;
    plot(db(fftshift(Tx_signal_pa)), 'DisplayName','NoDPD');
    hold on;
%     plot(db(fftshift(Oq_ruilc_sig_pa)), 'DisplayName','OQ RUILC');
%     plot(db(fftshift(Oq_ilc_sig_pa)), 'DisplayName','OQ ILC');
%     plot(db(fftshift(Ruilc_sig_pa)), 'DisplayName','RUILC');
%     plot(db(fftshift(Uilc_sig_pa)), 'DisplayName','UILC');
    plot(db(fftshift(Oqc_ilc_sig_pa)), 'DisplayName','OQC ILC');
    plot(db(fftshift(Ilc_sig_pa)), 'DisplayName','ILC');
%     plot(db(fftshift(Dl_sig_pa)), 'DisplayName','DL');
%     plot(db(fftshift(Idl_sig_pa)), 'DisplayName','IDL');
    hold off;
    legend('show');
    title('Magnitude spectra');
    xlabel('Frequency');
    ylabel('Magnitude (dB)');
    axis([-inf inf, -80 70]);
    
    if c_num_iter > 1
        % NMSE evolution
        % subplot('Position',[0.1, 0.05, 0.85, 0.25]);
        subplot(2,1,2);
        plot(nodpd_nmserr(1:iter), '-x', 'DisplayName','NoDPD');
        hold on;
        plot(oqc_ilc_nmserr(1:iter), '-x', 'DisplayName','OQC ILC');
        plot(oq_ruilc_nmserr(1:iter), '-x', 'DisplayName','OQ RUILC');
        plot(oq_ilc_nmserr(1:iter), '-x', 'DisplayName','OQ ILC');
    %     plot(ruilc_nmserr(1:iter), '-x', 'DisplayName','RUILC');
    %     plot(uilc_nmserr(1:iter), '-x', 'DisplayName','UILC');
        plot(ilc_nmserr(1:iter), '-x', 'DisplayName','ILC');
    %     plot(dl_nmserr(1:iter), '-x', 'DisplayName','DL');
    %     plot(idl_nmserr(1:iter), '-x', 'DisplayName','IDL');
        hold off;
        legend('show');
        title('Normalised mean square error');
        xlabel('Iteration cycle');
        ylabel('NMSE (dB)');
        axis([0 c_num_iter, -60 0]);
    else
        c_label_space = 5;
        y_pos = 20;
        text(0, y_pos, sprintf(' NoDPD NMSE = %.1f', nodpd_nmserr(1))); y_pos = y_pos - c_label_space;
        text(0, y_pos, sprintf(' OQC ILC NMSE = %.1f', oqc_ilc_nmserr(1))); y_pos = y_pos - c_label_space;
        text(0, y_pos, sprintf(' OQ RUILC NMSE = %.1f', oq_ruilc_nmserr(1))); y_pos = y_pos - c_label_space;
        text(0, y_pos, sprintf(' OQ ILC NMSE = %.1f', oq_ilc_nmserr(1))); y_pos = y_pos - c_label_space;
        text(0, y_pos, sprintf(' ILC NMSE = %.1f', ilc_nmserr(1))); y_pos = y_pos - c_label_space;
    end;
    
    pause(0.05);
end

% one-quadrature comparator iterative learning control
% comparator level is set from the low-speed DAC
% the digital circuit extracts the time when feedback signal crosses the
% set value, this value is afterwards given into DPD calculation 
function [tx_signal_pa, nmserr] = oqc_ilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M c_pa_gain c_pa_wp c_comp_val c_is_prob_sel;
    
    comp_val = c_comp_val;
    
    % generate low speed comparator voltages
    N = ceil(length(tx_signal)/length(comp_val));
    comp_sig = reshape(repmat(comp_val,N,1),length(comp_val)*N,1);
    comp_sig = comp_sig(1:length(tx_signal));        
    
    % calculate PA ouput
    tx_signal_pa = PA_Model((tx_signal*c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa))*max(abs(tx_signal));        
    %tx_signal_pa = abs(lscov(tx_signal_pa, tx_signal))*tx_signal_pa;
    tx_signal_pa = tx_signal_pa ./ c_pa_gain ./ c_pa_wp;        % normalize PA output to 1
    fb_sig = real(tx_signal_pa);
    
    % plot histogram of the feedback signal
    figure(2);
    subplot(2,2,1);
%     hist(fb_sig, linspace(-1/c_pa_gain,1/c_pa_gain,21));
%     title('Real feedback histogram');
%     axis([-1/c_pa_gain 1/c_pa_gain, 0, inf]);
    
    % plot IQ diagram here
    plot(tx_signal);
    hold on;
    plot(tx_signal_pa);
    for i = 1:length(comp_val)
        plot([comp_val(i) comp_val(i)], [-1 1], 'k');
    end
    hold off;
    axis([-2 2 -2 2]);
    title('IQ diagram');

    % compare the output of the PA with the comparison signal
    comp_out = real(fb_sig) > comp_sig;
    
    % detect edges on the comparator output
    comp_edges = [comp_out(2:end)-comp_out(1:end-1); 0];
    edge_ind = find(abs(comp_edges)==1);
    samples_oqc = length(edge_ind)
    
    if c_is_prob_sel
        % point preselection based on probability criteria
        % calculate expected absolute amplitude of the point based on the
        % amplitude of the input signal
        pt_amplitudes = abs(tx_signal(edge_ind));
        pt_rand = rand(size(pt_amplitudes));        % for every point from comparator, calculate random value
        pt_amp_lim = [0.2 0.35 0.5];
        pt_probabs_lim = [0.02 0.05 0.1 1];
        pt_probabs = zeros(size(pt_amplitudes));
        for i = 1:length(pt_amplitudes)
            % for every point select the probability at which the point is
            % selected into following calculations
            pt_probabs(i) = pt_probabs_lim(sum(pt_amp_lim < pt_amplitudes(i))+1);
        end
        edge_ind = edge_ind(pt_rand < pt_probabs);
        selected_samples_oqc = length(edge_ind)
    end
    
    % calculate exact edge time using linear interpolation
    edge_time = zeros(size(edge_ind));
    for i = 1:length(edge_ind)
        ei = edge_ind(i);
        edge_time(i) = ei + (fb_sig(ei)-comp_sig(ei))/(fb_sig(ei)-fb_sig(ei+1));
    end
    
    % plot the AM/AM and AM/PM characteristics with separated measured
    % points
    % interpolate the imaginary component to get precise value of the
    % imaginary value at the sample time
    fb_sig_i = imag(tx_signal_pa);
    fb_sig_dpd = zeros(size(edge_ind));
    tx_sig_dpd = zeros(size(edge_ind));
    for i = 1:length(edge_ind)
        ei = edge_ind(i);
        t = edge_time(i) - edge_ind(i);
        fb_sig_dpd(i) = comp_sig(ei) + 1j*(fb_sig_i(ei)+(fb_sig_i(ei+1)-fb_sig_i(ei))*t);
        tx_sig_dpd(i) = tx_signal(ei)+(tx_signal(ei+1)-tx_signal(ei))*t;
    end
    
    figure(2);
    subplot(2,2,3);
    plot(abs(tx_signal), angle(tx_signal)-angle(tx_signal_pa), '.',  'DisplayName','All samples');
    hold on;
    pm_diff = angle(tx_sig_dpd)-angle(fb_sig_dpd);
    pm_diff = pm_diff - 2*pi*(pm_diff > pi);
    pm_diff = pm_diff + 2*pi*(pm_diff < pi);
    plot(abs(tx_sig_dpd), pm_diff, '.',  'DisplayName','DPD');
    hold off;
    title('PM/AM characteristics');
    ylabel('Phase (rad)');
    legend('show', 'Location','southeast');
    
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_TimedMatrix(tx_signal, dpd_par.P, dpd_par.M, edge_time);
    M = [real(U) -imag(U)];
    x = comp_sig(edge_ind);
    b = ((M'*M)\(M'*x));
    pa_coefs = b(1:size(b,1)/2) + 1j*b(size(b,1)/2+1:end);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    %tx_signal_pa_mod = tx_signal_pa_mod./max(abs(tx_signal_pa_mod))*max(abs(tx_signal)); % normalize
    %tx_signal_pa_mod = abs(lscov(tx_signal_pa_mod, tx_signal))*tx_signal_pa_mod;
    %tx_signal_pa_mod = tx_signal_pa_mod / c_pa_wp;
    
    % plot AM/AM characteristisc
    figure(2);
    subplot(1,2,2);
    plot(abs(tx_signal), abs(tx_signal_pa_mod), '.',  'DisplayName','PA model');
    hold on;    
    plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','All samples');
    plot(abs(tx_sig_dpd), abs(fb_sig_dpd), '.',  'DisplayName','Selected samples');
    hold off;
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal);
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model((dpd_out * c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa))*max(abs(tx_signal)); % normalize
    %tx_signal_pa = abs(lscov(tx_signal_pa, tx_signal))*tx_signal_pa;
    tx_signal_pa = tx_signal_pa ./ c_pa_gain ./ c_pa_wp;
    
    % plot AM/AM characteristisc
    figure(2);
    subplot(1,2,2);
    hold on;    
    plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
    plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
    hold off;
    title('AM/AM characteristics');
    legend('show', 'Location','southeast');    
       
    % calculate NMSE
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    nmserr = nmse(tx_signal, tx_signal_pa);
end

% one-quadrature randomly undersampled iterative learning control
% randomly sampled (not at the specific sampling period rather with random
% sampling period - multiple of original sampling period) from feedback are taken for PA modeling
function [tx_signal_pa, nmserr] = oq_ruilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;
    
    c_undersampling_fact = 55;
    
    % generate random undersampling vector
    us = rand(size(tx_signal)) < (1/c_undersampling_fact);
    samples_oq = sum(us)
        
    % calculate PA ouput
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa));                         % normalize
    under_tx_signal_pa = tx_signal_pa(us);                                          % undersample
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_UndersampledMatrix(tx_signal./max(abs(tx_signal)), dpd_par.P, dpd_par.M, us);
    
    M = [real(U) -imag(U)];
    %x = real(rxSignal)./max(abs(real(rxSignal))); %normalize with respect to 1
    x = real(under_tx_signal_pa);
    b = ((M'*M)\(M'*x));
    pa_coefs = b(1:size(b,1)/2) + 1j*b(size(b,1)/2+1:end);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    tx_signal_pa_mod = tx_signal_pa_mod./max(abs(tx_signal_pa_mod)); % normalize with respect to 1
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal./max(abs(tx_signal)));
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end

% one-quadrature randomly undersampled iterative learning control
% randomly sampled (not at the specific sampling period rather with random
% sampling period - multiple of original sampling period) from feedback are taken for PA modeling
function [tx_signal_pa, nmserr] = oq_ilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;
    
    % calculate PA ouput
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa));                         % normalize
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_Matrix(tx_signal./max(abs(tx_signal)), dpd_par.P, dpd_par.M);
    
    M = [real(U) -imag(U)];
    %x = real(rxSignal)./max(abs(real(rxSignal))); %normalize with respect to 1
    x = real(tx_signal_pa);
    b = ((M'*M)\(M'*x));
    pa_coefs = b(1:size(b,1)/2) + 1j*b(size(b,1)/2+1:end);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    tx_signal_pa_mod = tx_signal_pa_mod./max(abs(tx_signal_pa_mod)); % normalize with respect to 1
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal./max(abs(tx_signal)));
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end

% randomly undersampled iterative learning control
% randomly sampled (not at the specific sampling period rather with random
% sampling period - multiple of original sampling period) from feedback are taken for PA modeling
function [tx_signal_pa, nmserr] = ruilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;
    
    c_undersampling_fact = 55;
    
    % generate random undersampling vector
    us = rand(size(tx_signal)) < (1/c_undersampling_fact);
        
    % calculate PA ouput
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    under_tx_signal_pa = tx_signal_pa(us);
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_UndersampledMatrix(tx_signal, dpd_par.P, dpd_par.M, us);
    pa_coefs = lscov(U, under_tx_signal_pa);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal);
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end


% undersampled iterative learning control
function [tx_signal_pa, nmserr] = uilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;
    
    c_undersampling_fact = 20;
    
    % calculate PA ouput
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa))*max(abs(tx_signal));        % normalize
    under_tx_signal_pa = tx_signal_pa(1:c_undersampling_fact:end);                  % undersample
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_UndersampledMatrix(tx_signal, dpd_par.P, dpd_par.M, c_undersampling_fact);
    pa_coefs = lscov(U, under_tx_signal_pa);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    tx_signal_pa_mod = tx_signal_pa_mod./max(abs(tx_signal_pa_mod))*max(abs(tx_signal)); % normalize
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal);
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end

% iterative learning control
function [tx_signal_pa, nmserr] = ilc_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M c_pa_wp c_pa_gain;
    
    % calculate PA output
    tx_signal_pa = PA_Model((tx_signal * c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = tx_signal_pa./max(abs(tx_signal_pa))*max(abs(tx_signal));    % normalize
    %tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;               % normalize
    tx_signal_pa = tx_signal_pa ./ c_pa_gain ./ c_pa_wp;
    
    % ILC: thus first calculate coeficients of the PA model
    U = DDR2_Matrix(tx_signal, dpd_par.P, dpd_par.M);
    pa_coefs = lscov(U, tx_signal_pa);
    
    % model PA and calculate its output
    tx_signal_pa_mod = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, pa_coefs).';
    
    % we have model output, calculate the required coefficients now
    U = DDR2_Matrix(tx_signal_pa_mod, dpd_par.P, dpd_par.M);
    coefs = lscov(U, tx_signal);
    
    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model((dpd_out*c_pa_wp).', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end

function [tx_signal_pa, nmserr, coefs] = dl_dpd(tx_signal, coefs, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;

    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
    
    % direct learning DPD coefficients - conventional DPD
    coefs = DirectLearning_DPD_DDR2(tx_signal, tx_signal_pa, coefs, dpd_par);
        
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end

function [tx_signal_pa, nmserr] = idl_dpd(tx_signal, dpd_par)
    global c_pa_coefs c_pa_P c_pa_M;
    
    % calculate DPD coefficients
    tx_signal_pa = PA_Model(tx_signal.', c_pa_coefs, c_pa_P, c_pa_M).';
    coefs = calc_DPD_DDR2(tx_signal_pa, tx_signal, dpd_par.P, dpd_par.M);

    % calculate DPD output
    dpd_out = DDR2_memory_polynomials(tx_signal, dpd_par.P, dpd_par.M, coefs).';

    % calculate PA model output
    tx_signal_pa = PA_Model(dpd_out.', c_pa_coefs, c_pa_P, c_pa_M).';
    %tx_signal_pa = (tx_signal_pa./max(abs(tx_signal_pa)));
    tx_signal_pa = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
       
    % calculate NMSE
    nmserr = nmse(tx_signal, tx_signal_pa);
end
