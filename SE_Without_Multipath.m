clc; clear; close all;

%% Common Parameters
N = 32;                 % Number of subcarriers/resources
M = 4;                  % QPSK modulation
k = log2(M);            % Bits per symbol
SNR_dB = 0:2:20;        % SNR range
numSymbols = 5000;      % Number of symbols per SNR point

% No multipath channel - set to flat channel
channelTaps = 1;        % Flat channel
maxDelay = 0;           % No delay
cp_len = 0;             % No cyclic prefix needed

% No equalization needed
equalizationType = 'None'; 
estimateChannel = 0;    % No channel estimation needed

%% Optical-Specific Parameters
% OFDM-IM Parameters (for DCO/ACO-OFDM-IM)
group_size = 8;         % Subcarriers per group
active_per_group = 4;   % Active subcarriers per group
numGroups = N / group_size;

% SCMA Parameters (for DCO/ACO-SCMA)
scma_users_optical = 2;         % Number of users
scma_resources_optical = 4;     % Number of resources
scma_codebooks_optical = {      % SCMA codebooks (real-valued for optical)
    [1, -1, 0.5, -0.5];  % User 1
    [0.5, -0.5, 1, -1]   % User 2
};

%% RF-Specific Parameters
% OFDM-IM Parameters for RF
group_size_rf = 4;         % Subcarriers per group
active_per_group_rf = 3;   % Active subcarriers per group
numGroups_rf = N / group_size_rf;

% SCMA Parameters for RF
scma_users_rf = 4;         % Number of users (150% overloading)
scma_resources_rf = 2;     % Number of resources
scma_codebooks_rf = cell(scma_users_rf, 1);
for u = 1:scma_users_rf
    scma_codebooks_rf{u} = exp(1j*2*pi*rand(M,1)); % Random codebooks
end

%% Preallocate Spectral Efficiency arrays
% Optical variants
se_DCOOFDM = zeros(size(SNR_dB));
se_ACOOFDM = zeros(size(SNR_dB));
se_DCOTFS = zeros(size(SNR_dB));
se_ACOTFS = zeros(size(SNR_dB));
se_DCOCDM = zeros(size(SNR_dB));
se_ACOCDM = zeros(size(SNR_dB));
se_DCOOFDMIM = zeros(size(SNR_dB));
se_ACOOFDMIM = zeros(size(SNR_dB));
se_DCOSCMA = zeros(size(SNR_dB));
se_ACOSCMA = zeros(size(SNR_dB));

% RF variants
se_OFDM = zeros(size(SNR_dB));
se_OTFS = zeros(size(SNR_dB));
se_OCDM = zeros(size(SNR_dB));
se_OFDMIM = zeros(size(SNR_dB));
se_SCMA = zeros(size(SNR_dB));

%% Main Simulation Loop
for snrIdx = 1:length(SNR_dB)
    snr = SNR_dB(snrIdx);
    fprintf('Processing SNR = %d dB...\n', snr);
    snr_lin = 10^(snr/10);
    
    % No channel estimation needed for flat channel
    H_est = ones(N, 1);  % Flat frequency response
    noiseVar = 10^(-snr/10);
    
    % Active subcarriers indices (for optical SE calculation)
    active_sc_DCO = 2:N/2; % For DCO schemes
    active_sc_ACO = 2:2:N/2; % For ACO schemes (only even indices)
    H_est_DCO = H_est(active_sc_DCO);
    H_est_ACO = H_est(active_sc_ACO);

    %% =================== Optical Variants ===================
    % DCO-OFDM
    dcBias = 3; % Approximate value
    tx_power = 1; % Normalized power
    dc_power = dcBias^2;
    effective_snr = snr_lin / (1 + dc_power/tx_power);
    capacity_per_sc = log2(1 + (abs(H_est_DCO).^2 * effective_snr));
    se_DCOOFDM(snrIdx) = mean(capacity_per_sc) * (N/2-1)/N;

    % ACO-OFDM
    tx_power = 1;
    effective_snr = snr_lin;
    capacity_per_sc = log2(1 + (abs(H_est_ACO).^2 * effective_snr));
    se_ACOOFDM(snrIdx) = mean(capacity_per_sc) * (N/4)/N;

    % DCO-OTFS
    dcBias = 3;
    tx_power = 1;
    dc_power = dcBias^2;
    effective_snr = snr_lin / (1 + dc_power/tx_power);
    H_est_OTFS = H_est_DCO; % No Doppler effect
    capacity_per_sc = log2(1 + (abs(H_est_OTFS).^2 * effective_snr));
    se_DCOTFS(snrIdx) = mean(capacity_per_sc) * (N/2-1)/N;

    % ACO-OTFS
    tx_power = 1;
    effective_snr = snr_lin;
    H_est_OTFS = H_est_ACO; % No Doppler effect
    capacity_per_sc = log2(1 + (abs(H_est_OTFS).^2 * effective_snr));
    se_ACOTFS(snrIdx) = mean(capacity_per_sc) * (N/4)/N;

    % DCO-OCDM
    dcBias = 3;
    tx_power = 1;
    dc_power = dcBias^2;
    effective_snr = snr_lin / (1 + dc_power/tx_power);
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N);
    chirp_response = abs(fft(chirp)).^2;
    chirp_response = chirp_response(active_sc_DCO);
    capacity_per_sc = log2(1 + (abs(H_est_DCO).^2 .* chirp_response * effective_snr));
    se_DCOCDM(snrIdx) = mean(capacity_per_sc) * (N/2-1)/N;

    % ACO-OCDM
    tx_power = 1;
    effective_snr = snr_lin;
    chirp_response = abs(fft(chirp)).^2;
    chirp_response = chirp_response(active_sc_ACO);
    capacity_per_sc = log2(1 + (abs(H_est_ACO).^2 .* chirp_response * effective_snr));
    se_ACOCDM(snrIdx) = mean(capacity_per_sc) * (N/4)/N;

    % DCO-OFDM-IM
    dcBias = 3;
    tx_power = 1;
    dc_power = dcBias^2;
    effective_snr = snr_lin / (1 + dc_power/tx_power);
    combinations = nchoosek(1:group_size, active_per_group);
    numCombinations = size(combinations, 1);
    bits_per_group = log2(numCombinations) + active_per_group * k;
    H_est_mag = abs(H_est_DCO).^2;
    effective_snr_IM = zeros(numGroups, 1);
    active_groups = floor(length(H_est_mag)/group_size);
    
    for g = 1:active_groups
        start_idx = (g-1)*group_size + 1;
        end_idx = g*group_size;
        group_indices = start_idx:min(end_idx, length(H_est_mag));
        
        if length(group_indices) == group_size
            active_snr = H_est_mag(group_indices) * effective_snr;
            effective_snr_IM(g) = mean(active_snr(combinations(1,:)));
        end
    end
    
    capacity_per_group = log2(1 + effective_snr_IM(1:active_groups));
    se_DCOOFDMIM(snrIdx) = mean(capacity_per_group) * (bits_per_group/group_size);

    % ACO-OFDM-IM
    tx_power = 1;
    effective_snr = snr_lin;
    H_est_mag = abs(H_est_ACO).^2;
    effective_snr_IM = zeros(numGroups, 1);
    active_groups = floor(length(H_est_mag)/group_size);
    
    for g = 1:active_groups
        start_idx = (g-1)*group_size + 1;
        end_idx = g*group_size;
        group_indices = start_idx:min(end_idx, length(H_est_mag));
        
        if length(group_indices) == group_size
            active_snr = H_est_mag(group_indices) * effective_snr;
            effective_snr_IM(g) = mean(active_snr(combinations(1,:)));
        end
    end
    
    capacity_per_group = log2(1 + effective_snr_IM(1:active_groups));
    se_ACOOFDMIM(snrIdx) = mean(capacity_per_group) * (bits_per_group/group_size);

    % DCO-SCMA
    dcBias = 3;
    tx_power = 1;
    dc_power = dcBias^2;
    effective_snr = snr_lin / (1 + dc_power/tx_power);
    H_est_scma = ones(scma_resources_optical, 1); % Flat channel
    capacity_per_resource = log2(1 + H_est_scma * effective_snr);
    se_DCOSCMA(snrIdx) = scma_users_optical * mean(capacity_per_resource);

    % ACO-SCMA
    tx_power = 1;
    effective_snr = snr_lin;
    H_est_scma = ones(scma_resources_optical, 1); % Flat channel
    capacity_per_resource = log2(1 + H_est_scma * effective_snr);
    se_ACOSCMA(snrIdx) = scma_users_optical * mean(capacity_per_resource);

    %% =================== RF Variants ===================
    % OFDM
    capacity_OFDM = mean(log2(1 + (abs(H_est).^2 * snr_lin)));
    se_OFDM(snrIdx) = capacity_OFDM;

    % OTFS
    effective_snr_otfs = (abs(H_est).^2)*snr_lin; % No Doppler effect
    se_OTFS(snrIdx) = mean(log2(1 + effective_snr_otfs));

    % OCDM
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N)/sqrt(N);
    chirp_response = abs(fft(chirp)).^2;
    effective_snr_ocdm = (abs(H_est).^2)*snr_lin.*chirp_response;
    se_OCDM(snrIdx) = mean(log2(1 + effective_snr_ocdm));

    % OFDM-IM
    combinations_rf = nchoosek(1:group_size_rf, active_per_group_rf);
    numCombinations_rf = size(combinations_rf, 1);
    if numCombinations_rf ~= 4 && numCombinations_rf ~= 8
        numCombinations_rf = 4;
        combinations_rf = combinations_rf(1:numCombinations_rf, :);
    end
    bits_per_group_rf = log2(numCombinations_rf) + active_per_group_rf * k;
    
    H_est_mag = abs(H_est).^2;
    effective_snr_im = zeros(numGroups_rf, 1);
    for g = 1:numGroups_rf
        idx = (g-1)*group_size_rf + (1:group_size_rf);
        active_snr = H_est_mag(idx) * snr_lin;
        effective_snr_im(g) = mean(active_snr(combinations_rf(1,:)));
    end
    capacity_per_group = log2(1 + effective_snr_im);
    se_OFDMIM(snrIdx) = mean(capacity_per_group) * (bits_per_group_rf/group_size_rf);

    % SCMA
    H_est_scma = ones(scma_resources_rf, 1); % Flat channel
    effective_snr_scma = H_est_scma * snr_lin;
    se_SCMA(snrIdx) = scma_users_rf * log2(1 + mean(effective_snr_scma));
end

%% Plotting
figure;
plot(SNR_dB, se_DCOOFDM, '-.o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM'); hold on;
plot(SNR_dB, se_DCOTFS, '-.x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
plot(SNR_dB, se_DCOCDM, '-.+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
plot(SNR_dB, se_DCOOFDMIM, '-.s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
plot(SNR_dB, se_DCOSCMA, '-.d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

% ACO variants (dashed lines)
plot(SNR_dB, se_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
plot(SNR_dB, se_ACOTFS, '--x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
plot(SNR_dB, se_ACOCDM, '--+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
plot(SNR_dB, se_ACOOFDMIM, '--s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
plot(SNR_dB, se_ACOSCMA, '--d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('SNR (dB)'); 
ylabel('Spectral Efficiency (bit/s/Hz)');
legend('Location', 'northwest', 'NumColumns', 2);
grid on;
set(gca, 'FontSize', 12);

% RF Variants
figure;
plot(SNR_dB, se_OFDM, '-o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'OFDM'); hold on;
plot(SNR_dB, se_OTFS, '-x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'OTFS');
plot(SNR_dB, se_OCDM, '-+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'OCDM');
plot(SNR_dB, se_OFDMIM, '-s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'OFDM-IM');
plot(SNR_dB, se_SCMA, '-d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'SCMA');

plot(SNR_dB, se_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
plot(SNR_dB, se_ACOTFS, '--x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
plot(SNR_dB, se_DCOCDM, '-.+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
plot(SNR_dB, se_ACOOFDMIM, '--s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
plot(SNR_dB, se_ACOSCMA, '--d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');
xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bit/s/Hz)');
legend('Location', 'northwest', 'NumColumns', 2);
grid on;
set(gca, 'FontSize', 12);