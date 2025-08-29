clc; clear; close all;

%% Parameters
N = 32;                 % Number of subcarriers/resources
M = 4;                  % QPSK modulation
k = log2(M);            % Bits per symbol
SNR_dB = 0:2:20;        % SNR range
numSymbols = 5000;      % Number of symbols per SNR point

% Multipath channel parameters
channelTaps = [0.8, 0.05, 0.1]; % Channel tap gains
maxDelay = length(channelTaps)-1; % Maximum delay in samples
cp_len = maxDelay;      % Cyclic prefix length

% Equalization parameters
equalizationType = 'MMSE'; % 'ZF' or 'MMSE'
estimateChannel = 1;     % Channel estimation flag

%% OFDM-IM Parameters
group_size = 8;         % Subcarriers per group
active_per_group = 4;   % Active subcarriers per group
numGroups = N / group_size;

% Generate all possible active subcarrier combinations
combinations = nchoosek(1:group_size, active_per_group);
numCombinations = size(combinations, 1);
if numCombinations ~= 4 && numCombinations ~= 8
    numCombinations = 4;
    combinations = combinations(1:numCombinations, :);
end
bits_per_group = log2(numCombinations) + active_per_group * k;

%% SCMA Parameters
scma_users = 2;         % Number of users
scma_resources = 4;     % Number of resources
scma_codebooks = {      % SCMA codebooks
    [1+1i, -1+1i, 1-1i, -1-1i]/sqrt(2);  % User 1
    [1+1i, 1-1i, -1+1i, -1-1i]/sqrt(2)   % User 2
};

%% Preallocate BER arrays
ber_OFDM = zeros(size(SNR_dB));
ber_DCOOFDM = zeros(size(SNR_dB));
ber_ACOOFDM = zeros(size(SNR_dB));
ber_OTFS = zeros(size(SNR_dB));
ber_DCOTFS = zeros(size(SNR_dB));
ber_ACOTFS = zeros(size(SNR_dB));
ber_OCDM = zeros(size(SNR_dB));
ber_OFDMIM = zeros(size(SNR_dB));
ber_SCMA = zeros(size(SNR_dB));

ber_DCOCDM = zeros(size(SNR_dB));
ber_ACOCDM = zeros(size(SNR_dB));

ber_DCOOFDMIM = zeros(size(SNR_dB));
ber_ACOOFDMIM = zeros(size(SNR_dB));

ber_DCOSCMA= zeros(size(SNR_dB));
ber_ACOSCMA = zeros(size(SNR_dB));
%% Main Simulation Loop
for snrIdx = 1:length(SNR_dB)
    snr = SNR_dB(snrIdx);
    fprintf('Processing SNR = %d dB...\n', snr);
    
    % Channel estimation
    if estimateChannel
        pilotSymbols = qammod(randi([0 M-1], N, 1), M, 'UnitAveragePower', true);
        pilotOFDM = ifft(pilotSymbols);
        pilotOFDM_CP = [pilotOFDM(end-cp_len+1:end); pilotOFDM];
        rxPilot = filter(channelTaps, 1, pilotOFDM_CP);
        rxPilot = rxPilot(cp_len+1:end);
        H_est = fft(rxPilot) ./ fft(pilotOFDM);
    else
        H_est = fft(channelTaps, N).';
    end
    noiseVar = 10^(-snr/10);

    %% =================== OFDM ===================
    bits_OFDM = randi([0 1], N*numSymbols*k, 1);
    symbols_OFDM = qammod(bits2int(bits_OFDM, k), M, 'UnitAveragePower', true);
    tx_OFDM = ifft(reshape(symbols_OFDM, N, numSymbols));
    
    % Add CP and transmit through channel
    tx_OFDM_CP = [tx_OFDM(end-cp_len+1:end, :); tx_OFDM];
    rx_OFDM_CP = zeros(size(tx_OFDM_CP));
    for sym = 1:numSymbols
       rx_OFDM_CP(:, sym) = filter(channelTaps, 1, tx_OFDM_CP(:, sym));
    end
    rx_OFDM_CP = awgn(rx_OFDM_CP, snr, 'measured');
    
    % Remove CP and equalize
    rx_OFDM = rx_OFDM_CP(cp_len+1:end, :);
    rxSymbols_OFDM = fft(rx_OFDM);
    if strcmp(equalizationType, 'ZF')
        rxSymbols_OFDM = rxSymbols_OFDM ./ H_est;
    else % MMSE
        rxSymbols_OFDM = rxSymbols_OFDM .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    end
    
    % Demodulate
    rxBits_OFDM = int2bits(qamdemod(rxSymbols_OFDM(:), M, 'UnitAveragePower', true), k);
    ber_OFDM(snrIdx) = mean(bits_OFDM ~= rxBits_OFDM(1:length(bits_OFDM)));

    %% =================== DCO-OFDM ===================
    bits_DCOOFDM = randi([0 1], (N/2-1)*numSymbols*k, 1);
    symbols_DCOOFDM = qammod(bits2int(bits_DCOOFDM, k), M, 'UnitAveragePower', true);
    data = reshape(symbols_DCOOFDM, (N/2-1), numSymbols);
    data2 = [zeros(1,numSymbols);
            data;
            zeros(1,numSymbols);
            flipud(conj(data))];
    tx_DCOOFDM = ifft(data2);
    dcBias = -min(real(tx_DCOOFDM(:)));
    tx_DCOOFDM = real(tx_DCOOFDM) + dcBias;

    % Add CP and transmit through channel
    tx_DCOOFDM_CP = [tx_DCOOFDM(end-cp_len+1:end, :); tx_DCOOFDM];
    rx_DCOOFDM_CP = zeros(size(tx_DCOOFDM_CP));
    for sym = 1:numSymbols
       rx_DCOOFDM_CP(:, sym) = filter(channelTaps, 1, tx_DCOOFDM_CP(:, sym));
    end
    rx_DCOOFDM_CP = awgn(rx_DCOOFDM_CP, snr, 'measured');
    
    % Remove CP and equalize
    rx_DCOOFDM = rx_DCOOFDM_CP(cp_len+1:end, :);
    rxSymbols_DCOOFDM = fft(rx_DCOOFDM);
    if strcmp(equalizationType, 'ZF')
        rxSymbols_DCOOFDM = rxSymbols_DCOOFDM ./ H_est;
    else % MMSE
        rxSymbols_DCOOFDM = rxSymbols_DCOOFDM .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    end
    
    % Demodulate
    rxBits_DCOOFDM = int2bits(qamdemod(rxSymbols_DCOOFDM(2:N/2,:), M, 'UnitAveragePower', true), k);
    ber_DCOOFDM(snrIdx) = mean(bits_DCOOFDM ~= rxBits_DCOOFDM(1:length(bits_DCOOFDM)));

   %% =================== ACO-OFDM ===================
    bits_ACOOFDM = randi([0 1], (N/4)*numSymbols*k, 1);
    symbols_ACOOFDM = qammod(bits2int(bits_ACOOFDM, k), M, 'UnitAveragePower', true);
    data = reshape(symbols_ACOOFDM, (N/4), numSymbols);
    data2 = zeros(N,numSymbols);
    for kk=1:N/2
        if ~rem(kk,2)
            data2(kk,:) = data(kk/2,:);
        end
    end
    data2(N/2+2:end,:) = flipud(conj(data2(2:N/2,:)));
    tx_ACOOFDM = ifft(data2);
    tx_ACOOFDM(tx_ACOOFDM<0) = 0;    %Clipping
    
    % Add CP and transmit through channel
    tx_ACOOFDM_CP = [tx_ACOOFDM(end-cp_len+1:end, :); tx_ACOOFDM];
    rx_ACOOFDM_CP = zeros(size(tx_ACOOFDM_CP));
    for sym = 1:numSymbols
       rx_ACOOFDM_CP(:, sym) = filter(channelTaps, 1, tx_ACOOFDM_CP(:, sym));
    end
    rx_ACOOFDM_CP = awgn(rx_ACOOFDM_CP, snr, 'measured');
    
    % Remove CP and equalize
    rx_ACOOFDM = rx_ACOOFDM_CP(cp_len+1:end, :);
    rxSymbols_ACOOFDM = fft(rx_ACOOFDM);
    if strcmp(equalizationType, 'ZF')
        rxSymbols_ACOOFDM = rxSymbols_ACOOFDM ./ H_est;
    else % MMSE
        rxSymbols_ACOOFDM = rxSymbols_ACOOFDM .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    end
    
    % Demodulate
    rxBits_ACOOFDM = int2bits(qamdemod(rxSymbols_ACOOFDM(2:2:N/2,:), M, 'UnitAveragePower', true), k);
    ber_ACOOFDM(snrIdx) = mean(bits_ACOOFDM ~= rxBits_ACOOFDM(1:length(bits_ACOOFDM)));

    %% =================== OTFS ===================
    bits_OTFS = randi([0 1], N*numSymbols*k, 1);
    symbols_OTFS = qammod(bits2int(bits_OTFS, k), M, 'UnitAveragePower', true);

    % Reshape symbols into delay-Doppler grid
    dd_grid = reshape(symbols_OTFS, N, numSymbols);

    % OTFS modulation (ISFFT)
    tx_OTFS = ifft(fft(dd_grid, [], 2), [], 1);

    % Add CP and transmit
    tx_OTFS_CP = [tx_OTFS(end-cp_len+1:end, :); tx_OTFS];
    rx_OTFS_CP = zeros(size(tx_OTFS_CP));
    for sym = 1:numSymbols
        rx_OTFS_CP(:, sym) = filter(channelTaps, 1, tx_OTFS_CP(:, sym));
    end
    rx_OTFS_CP = awgn(rx_OTFS_CP, snr, 'measured');

    % Remove CP and OTFS demodulation (SFFT)
    rx_OTFS = rx_OTFS_CP(cp_len+1:end, :);
    rxSymbols_OTFS = ifft(fft(rx_OTFS, [], 1), [], 2);

    % Demodulate
    rxBits_OTFS = int2bits(qamdemod(rxSymbols_OTFS(:), M, 'UnitAveragePower', true), k);
    ber_OTFS(snrIdx) = mean(bits_OTFS ~= rxBits_OTFS(1:length(bits_OTFS)));

%% =================== DCO-OTFS ===================
bits_DCOTFS = randi([0 1], (N/2-1)*numSymbols*k, 1);
symbols_DCOTFS = qammod(bits2int(bits_DCOTFS, k), M, 'UnitAveragePower', true);

% Create Hermitian symmetric input for real output
data = reshape(symbols_DCOTFS, (N/2-1), numSymbols);
data2 = [zeros(1,numSymbols);  % DC subcarrier (must be zero)
         data; 
         zeros(1,numSymbols);  % Nyquist frequency
         flipud(conj(data))];  % Hermitian symmetry

% Reshape into delay-Doppler grid with proper structure
dd_grid = zeros(N, numSymbols);
dd_grid(2:N/2,:) = data2(2:N/2,:);          % Positive frequencies
dd_grid(N/2+2:end,:) = conj(flipud(data2(2:N/2,:))); % Negative frequencies

% OTFS modulation (ISFFT)
tx_DCOTFS = ifft(fft(dd_grid, [], 2), [], 1, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_DCOTFS), 'Transmitted signal is not real-valued');

% Add DC bias to make signal non-negative
dcBias = -min(tx_DCOTFS(:));
tx_DCOTFS = tx_DCOTFS + dcBias;

% Add CP and transmit through channel
tx_DCOTFS_CP = [tx_DCOTFS(end-cp_len+1:end, :); tx_DCOTFS]
rx_DCOTFS_CP = zeros(size(tx_DCOTFS_CP));
for sym = 1:numSymbols
    rx_DCOTFS_CP(:, sym) = filter(channelTaps, 1, tx_DCOTFS_CP(:, sym));
end
rx_DCOTFS_CP = awgn(rx_DCOTFS_CP, snr, 'measured');

% Remove CP and OTFS demodulation (SFFT)
rx_DCOTFS = rx_DCOTFS_CP(cp_len+1:end, :);
rxSymbols_DCOTFS = ifft(fft(rx_DCOTFS, [], 1), [], 2);

% Remove DC bias effect and extract data
rxSymbols_DCOTFS = rxSymbols_DCOTFS(2:N/2,:);

% Demodulate
rxBits_DCOTFS = int2bits(qamdemod(rxSymbols_DCOTFS(:), M, 'UnitAveragePower', true), k);
ber_DCOTFS(snrIdx) = mean(bits_DCOTFS ~= rxBits_DCOTFS(1:length(bits_DCOTFS)));

%% =================== ACO-OTFS ===================
bits_ACOTFS = randi([0 1], (N/4)*numSymbols*k, 1);
symbols_ACOTFS = qammod(bits2int(bits_ACOTFS, k), M, 'UnitAveragePower', true);

% Structure input for ACO-OTFS (only odd subcarriers carry information)
data = reshape(symbols_ACOTFS, (N/4), numSymbols);
dd_grid = zeros(N, numSymbols);
for kk = 1:N/2
    if mod(kk,2) == 0  % Only even indices (which become odd after adding DC)
        dd_grid(kk,:) = data(kk/2,:);
    end
end
dd_grid(N/2+2:end,:) = conj(flipud(dd_grid(2:N/2,:))); % Hermitian symmetry

% OTFS modulation (ISFFT)
tx_ACOTFS = ifft(fft(dd_grid, [], 2), [], 1, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_ACOTFS), 'Transmitted signal is not real-valued');

% Clip negative parts to zero
tx_ACOTFS(tx_ACOTFS < 0) = 0;

% Add CP and transmit
tx_ACOTFS_CP = [tx_ACOTFS(end-cp_len+1:end, :); tx_ACOTFS];
rx_ACOTFS_CP = zeros(size(tx_ACOTFS_CP));
for sym = 1:numSymbols
    rx_ACOTFS_CP(:, sym) = filter(channelTaps, 1, tx_ACOTFS_CP(:, sym));
end
rx_ACOTFS_CP = awgn(rx_ACOTFS_CP, snr, 'measured');

% Remove CP and OTFS demodulation (SFFT)
rx_ACOTFS = rx_ACOTFS_CP(cp_len+1:end, :);
rxSymbols_ACOTFS = ifft(fft(rx_ACOTFS, [], 1), [], 2);

% Extract data from even subcarriers (which carried the information)
rxSymbols_ACOTFS = rxSymbols_ACOTFS(2:2:N/2,:);

% Demodulate
rxBits_ACOTFS = int2bits(qamdemod(rxSymbols_ACOTFS(:), M, 'UnitAveragePower', true), k);
ber_ACOTFS(snrIdx) = mean(bits_ACOTFS ~= rxBits_ACOTFS(1:length(bits_ACOTFS)));
    %% =================== OCDM ===================
    bits_OCDM = randi([0 1], N*numSymbols*k, 1);
    symbols_OCDM = qammod(bits2int(bits_OCDM, k), M, 'UnitAveragePower', true);
    
    % OCDM modulation (DAFT)
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N);
    tx_OCDM = ifft(reshape(symbols_OCDM, N, numSymbols) .* chirp);
    
    % Add CP and transmit
    tx_OCDM_CP = [tx_OCDM(end-cp_len+1:end, :); tx_OCDM];
    rx_OCDM_CP = zeros(size(tx_OCDM_CP));
    for sym = 1:numSymbols
        rx_OCDM_CP(:, sym) = filter(channelTaps, 1, tx_OCDM_CP(:, sym));
    end
    rx_OCDM_CP = awgn(rx_OCDM_CP, snr, 'measured');
    
    % Remove CP and demodulate
    rx_OCDM = rx_OCDM_CP(cp_len+1:end, :);
    y = fft(rx_OCDM) ./ chirp;
    
    % Equalization
    if strcmp(equalizationType, 'ZF')
        y = y ./ H_est;
    else % MMSE
        y = y .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    end
    
    % Demodulate
    rxBits_OCDM = int2bits(qamdemod(y(:), M, 'UnitAveragePower', true), k);
    ber_OCDM(snrIdx) = mean(bits_OCDM ~= rxBits_OCDM(1:length(bits_OCDM)));

    %% =================== DCO-OCDM ===================
bits_DCOCDM = randi([0 1], (N/2-1)*numSymbols*k, 1);
symbols_DCOCDM = qammod(bits2int(bits_DCOCDM, k), M, 'UnitAveragePower', true);

% Create Hermitian symmetric input for real output
data = reshape(symbols_DCOCDM, (N/2-1), numSymbols);
data2 = [zeros(1,numSymbols);  % DC subcarrier (must be zero)
         data; 
         zeros(1,numSymbols);   % Nyquist frequency
         flipud(conj(data))];  % Hermitian symmetry

% Reshape into proper structure for OCDM
input_sym = zeros(N, numSymbols);
input_sym(2:N/2,:) = data2(2:N/2,:);          % Positive frequencies
input_sym(N/2+2:end,:) = conj(flipud(data2(2:N/2,:))); % Negative frequencies

% OCDM modulation (DAFT)
n = (0:N-1).';
chirp = exp(1j*pi*(n.^2)/N);
tx_DCOCDM = ifft(input_sym .* chirp, [], 1, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_DCOCDM), 'Transmitted signal is not real-valued');

% Add DC bias to make signal non-negative
dcBias = -min(tx_DCOCDM(:));
tx_DCOCDM = tx_DCOCDM + dcBias;

% Add CP and transmit through channel
tx_DCOCDM_CP = [tx_DCOCDM(end-cp_len+1:end, :); tx_DCOCDM];
rx_DCOCDM_CP = zeros(size(tx_DCOCDM_CP));
for sym = 1:numSymbols
    rx_DCOCDM_CP(:, sym) = filter(channelTaps, 1, tx_DCOCDM_CP(:, sym));
end
rx_DCOCDM_CP = awgn(rx_DCOCDM_CP, snr, 'measured');

% Remove CP and demodulate
rx_DCOCDM = rx_DCOCDM_CP(cp_len+1:end, :);
y = fft(rx_DCOCDM) ./ chirp;

% Equalization
if strcmp(equalizationType, 'ZF')
    y = y ./ H_est;
else % MMSE
    y = y .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
end

% Remove DC bias effect and extract data
y = y(2:N/2,:);

% Demodulate
rxBits_DCOCDM = int2bits(qamdemod(y(:), M, 'UnitAveragePower', true), k);
ber_DCOCDM(snrIdx) = mean(bits_DCOCDM ~= rxBits_DCOCDM(1:length(bits_DCOCDM)));

%% =================== ACO-OCDM ===================
bits_ACOCDM = randi([0 1], (N/4)*numSymbols*k, 1);
symbols_ACOCDM = qammod(bits2int(bits_ACOCDM, k), M, 'UnitAveragePower', true);

% Structure input for ACO-OCDM (only odd subcarriers carry information)
data = reshape(symbols_ACOCDM, (N/4), numSymbols);
input_sym = zeros(N, numSymbols);
for kk = 1:N/2
    if mod(kk,2) == 0  % Only even indices (which become odd after adding DC)
        input_sym(kk,:) = data(kk/2,:);
    end
end
input_sym(N/2+2:end,:) = conj(flipud(input_sym(2:N/2,:))); % Hermitian symmetry

% OCDM modulation (DAFT)
n = (0:N-1).';
chirp = exp(1j*pi*(n.^2)/N);
tx_ACOCDM = ifft(input_sym .* chirp, [], 1, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_ACOCDM), 'Transmitted signal is not real-valued');

% Clip negative parts to zero
tx_ACOCDM(tx_ACOCDM < 0) = 0;

% Add CP and transmit
tx_ACOCDM_CP = [tx_ACOCDM(end-cp_len+1:end, :); tx_ACOCDM];
rx_ACOCDM_CP = zeros(size(tx_ACOCDM_CP));
for sym = 1:numSymbols
    rx_ACOCDM_CP(:, sym) = filter(channelTaps, 1, tx_ACOCDM_CP(:, sym));
end
rx_ACOCDM_CP = awgn(rx_ACOCDM_CP, snr, 'measured');

% Remove CP and demodulate
rx_ACOCDM = rx_ACOCDM_CP(cp_len+1:end, :);
y = fft(rx_ACOCDM) ./ chirp;

% Equalization
if strcmp(equalizationType, 'ZF')
    y = y ./ H_est;
else % MMSE
    y = y .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
end

% Extract data from even subcarriers (which carried the information)
y = y(2:2:N/2,:);

% Demodulate
rxBits_ACOCDM = int2bits(qamdemod(y(:), M, 'UnitAveragePower', true), k);
ber_ACOCDM(snrIdx) = mean(bits_ACOCDM ~= rxBits_ACOCDM(1:length(bits_ACOCDM)));

    %% =================== OFDM-IM ===================
    totalBits_OFDMIM = numGroups*bits_per_group*numSymbols;
    bits_OFDMIM = randi([0 1], totalBits_OFDMIM, 1);
    
    % Reshape bits with padding if needed
    padBits = mod(length(bits_OFDMIM), bits_per_group);
    if padBits > 0
        bits_OFDMIM = [bits_OFDMIM; zeros(bits_per_group-padBits, 1)];
    end
    
    bits_reshaped = reshape(bits_OFDMIM, bits_per_group, []).';
    index_bits = bits_reshaped(:, 1:log2(numCombinations));
    data_bits = bits_reshaped(:, log2(numCombinations)+1:end);
    index_decimal = bi2de(index_bits, 'left-msb') + 1;

    % Convert data bits to integers
    data_ints = zeros(size(data_bits,1), active_per_group);
    for i = 1:active_per_group
        start_bit = (i-1)*k + 1;
        end_bit = i*k;
        data_ints(:,i) = bits2int(data_bits(:, start_bit:end_bit), k);
    end
    
    % Modulation
    mod_data = qammod(data_ints(:), M, 'UnitAveragePower', true);
    mod_data = reshape(mod_data, [], active_per_group);

    % Create transmit blocks
    tx_blocks = zeros(N, numSymbols);
    for sym = 1:numSymbols
        for g = 1:numGroups
            idx = (g-1)*group_size + combinations(index_decimal((sym-1)*numGroups + g), :);
            tx_blocks(idx, sym) = mod_data((sym-1)*numGroups + g, :);
        end
    end

    % IFFT and add CP
    tx_OFDMIM = ifft(tx_blocks, N);
    tx_OFDMIM_CP = [tx_OFDMIM(end-cp_len+1:end, :); tx_OFDMIM];
    
    % Channel and noise
    rx_OFDMIM_CP = zeros(size(tx_OFDMIM_CP));
    for sym = 1:numSymbols
        rx_OFDMIM_CP(:, sym) = filter(channelTaps, 1, tx_OFDMIM_CP(:, sym));
    end
    rx_OFDMIM_CP = awgn(rx_OFDMIM_CP, snr, 'measured');
    
    % Remove CP and FFT
    rx_OFDMIM = rx_OFDMIM_CP(cp_len+1:end, :);
    rx_blocks = fft(rx_OFDMIM, N);
    
    % MMSE Equalization
    rx_blocks = rx_blocks .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);

    % OFDM-IM detection
    rx_bits_IM = [];
    for sym = 1:numSymbols
        for g = 1:numGroups
            idxs = (g-1)*group_size + (1:group_size);
            rx_group = rx_blocks(idxs, sym);
            
            % Find active subcarriers
            [~, sorted_idx] = sort(abs(rx_group).^2, 'descend');
            detected_active = sort(sorted_idx(1:active_per_group))';
            
            % Find matching combination
            diff = abs(combinations - detected_active);
            [~, comb_idx] = min(sum(diff, 2));
            
            % Decode index bits
            rx_index_bits = de2bi(comb_idx-1, log2(numCombinations), 'left-msb');
            
            % Decode data symbols
            rx_symbols = rx_group(combinations(comb_idx, :));
            rx_data_bits = int2bits(qamdemod(rx_symbols, M, 'UnitAveragePower', true), k);
            
            rx_bits_IM = [rx_bits_IM; rx_index_bits(:); rx_data_bits(:)];
        end
    end
    
    % Remove padding if any
    if padBits > 0
        rx_bits_IM = rx_bits_IM(1:end-(bits_per_group-padBits));
    end
    
    ber_OFDMIM(snrIdx) = mean(bits_OFDMIM ~= rx_bits_IM(1:length(bits_OFDMIM)));
    %% =================== DCO-OFDM-IM ===================

effective_N = N/2-1;  % Number of usable subcarriers (excluding DC and Nyquist)
effective_numGroups = floor(effective_N/group_size);
totalBits_DCOOFDMIM = effective_numGroups*bits_per_group*numSymbols;
bits_DCOOFDMIM = randi([0 1], totalBits_DCOOFDMIM, 1);

% Reshape bits with padding if needed
padBits = mod(length(bits_DCOOFDMIM), bits_per_group);
if padBits > 0
    bits_DCOOFDMIM = [bits_DCOOFDMIM; zeros(bits_per_group-padBits, 1)];
end

bits_reshaped = reshape(bits_DCOOFDMIM, bits_per_group, []).';
index_bits = bits_reshaped(:, 1:log2(numCombinations));
data_bits = bits_reshaped(:, log2(numCombinations)+1:end);
index_decimal = bi2de(index_bits, 'left-msb') + 1;

% Convert data bits to integers
data_ints = zeros(size(data_bits,1), active_per_group);
for i = 1:active_per_group
    start_bit = (i-1)*k + 1;
    end_bit = i*k;
    data_ints(:,i) = bits2int(data_bits(:, start_bit:end_bit), k);
end

% Modulation
mod_data = qammod(data_ints(:), M, 'UnitAveragePower', true);
mod_data = reshape(mod_data, [], active_per_group);

% Create transmit blocks with Hermitian symmetry
tx_blocks = zeros(N, numSymbols);
for sym = 1:numSymbols
    for g = 1:effective_numGroups
        % Map to positive frequencies (skip DC at index 1)
        idx = (g-1)*group_size + combinations(index_decimal((sym-1)*effective_numGroups + g), :) + 1;
        tx_blocks(idx, sym) = mod_data((sym-1)*effective_numGroups + g, :);
    end
    % Create Hermitian symmetric counterpart
    tx_blocks(N/2+2:end, sym) = conj(flipud(tx_blocks(2:N/2, sym)));
end

% OFDM-IM modulation
tx_DCOOFDMIM = ifft(tx_blocks, N, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_DCOOFDMIM), 'Transmitted signal is not real-valued');

% Add DC bias to make signal non-negative
dcBias = -min(tx_DCOOFDMIM(:));
tx_DCOOFDMIM = tx_DCOOFDMIM + dcBias;

% Add CP and transmit through channel
tx_DCOOFDMIM_CP = [tx_DCOOFDMIM(end-cp_len+1:end, :); tx_DCOOFDMIM]
rx_DCOOFDMIM_CP = zeros(size(tx_DCOOFDMIM_CP));
for sym = 1:numSymbols
    rx_DCOOFDMIM_CP(:, sym) = filter(channelTaps, 1, tx_DCOOFDMIM_CP(:, sym));
end
rx_DCOOFDMIM_CP = awgn(rx_DCOOFDMIM_CP, snr, 'measured');

% Remove CP and FFT
rx_DCOOFDMIM = rx_DCOOFDMIM_CP(cp_len+1:end, :);
rx_blocks = fft(rx_DCOOFDMIM);

% MMSE Equalization
rx_blocks = rx_blocks .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);

% DCO-OFDM-IM detection
rx_bits_DCOIM = [];
for sym = 1:numSymbols
    for g = 1:effective_numGroups
        idxs = (g-1)*group_size + (1:group_size) + 1; % +1 to skip DC
        rx_group = rx_blocks(idxs, sym);
        
        % Find active subcarriers
        [~, sorted_idx] = sort(abs(rx_group).^2, 'descend');
        detected_active = sort(sorted_idx(1:active_per_group))';
        
        % Find matching combination
        diff = abs(combinations - detected_active);
        [~, comb_idx] = min(sum(diff, 2));
        
        % Decode index bits
        rx_index_bits = de2bi(comb_idx-1, log2(numCombinations), 'left-msb');
        
        % Decode data symbols
        rx_symbols = rx_group(combinations(comb_idx, :));
        rx_data_bits = int2bits(qamdemod(rx_symbols, M, 'UnitAveragePower', true), k);
        
        rx_bits_DCOIM = [rx_bits_DCOIM; rx_index_bits(:); rx_data_bits(:)];
    end
end

% Remove padding if any
if padBits > 0
    rx_bits_DCOIM = rx_bits_DCOIM(1:end-(bits_per_group-padBits));
end

ber_DCOOFDMIM(snrIdx) = mean(bits_DCOOFDMIM ~= rx_bits_DCOIM(1:length(bits_DCOOFDMIM)));

  %% =================== ACO-OFDM-IM ===================
numEvenSubcarriers = N/2;  % Number of even-indexed subcarriers
effectiveGroupsACO = floor(numEvenSubcarriers/group_size); % Number of complete groups
bitsPerGroupACO = log2(numCombinations) + active_per_group*k;
totalBits_ACOOFDMIM = effectiveGroupsACO * bitsPerGroupACO * numSymbols;

% Ensure totalBits_ACOOFDMIM is integer
totalBits_ACOOFDMIM = floor(totalBits_ACOOFDMIM);
bits_ACOOFDMIM = randi([0 1], totalBits_ACOOFDMIM, 1);

% Reshape bits with padding if needed
padBits = mod(length(bits_ACOOFDMIM), bitsPerGroupACO);
if padBits > 0
    bits_ACOOFDMIM = [bits_ACOOFDMIM; zeros(bitsPerGroupACO-padBits, 1)];
end

bits_reshaped = reshape(bits_ACOOFDMIM, bitsPerGroupACO, []).';
index_bits = bits_reshaped(:, 1:log2(numCombinations));
data_bits = bits_reshaped(:, log2(numCombinations)+1:end);
index_decimal = bi2de(index_bits, 'left-msb') + 1;

% Convert data bits to integers
data_ints = zeros(size(data_bits,1), active_per_group);
for i = 1:active_per_group
    start_bit = (i-1)*k + 1;
    end_bit = i*k;
    data_ints(:,i) = bits2int(data_bits(:, start_bit:end_bit), k);
end

% Modulation
mod_data = qammod(data_ints(:), M, 'UnitAveragePower', true);
mod_data = reshape(mod_data, [], active_per_group);

% Create transmit blocks using only even subcarriers (2,4,6,...)
tx_blocks = zeros(N, numSymbols);
for sym = 1:numSymbols
    for g = 1:effectiveGroupsACO
        % Calculate even subcarrier indices for this group
        startSubcarrier = 2 + (g-1)*group_size;
        activeSubcarriers = startSubcarrier + combinations(index_decimal((sym-1)*effectiveGroupsACO + g), :) - 1;
        
        % Map data to active subcarriers
        tx_blocks(activeSubcarriers, sym) = mod_data((sym-1)*effectiveGroupsACO + g, :);
    end
    % Create Hermitian symmetric counterpart
    tx_blocks(N/2+2:end, sym) = conj(flipud(tx_blocks(2:N/2, sym)));
end

% OFDM-IM modulation
tx_ACOOFDMIM = ifft(tx_blocks, N, 'symmetric');

% Verify real-valued signal
assert(isreal(tx_ACOOFDMIM), 'Transmitted signal is not real-valued');

% Clip negative parts to zero
tx_ACOOFDMIM(tx_ACOOFDMIM < 0) = 0;

% Add CP and transmit
tx_ACOOFDMIM_CP = [tx_ACOOFDMIM(end-cp_len+1:end, :); tx_ACOOFDMIM];
rx_ACOOFDMIM_CP = zeros(size(tx_ACOOFDMIM_CP));
for sym = 1:numSymbols
    rx_ACOOFDMIM_CP(:, sym) = filter(channelTaps, 1, tx_ACOOFDMIM_CP(:, sym));
end
rx_ACOOFDMIM_CP = awgn(rx_ACOOFDMIM_CP, snr, 'measured');

% Remove CP and FFT
rx_ACOOFDMIM = rx_ACOOFDMIM_CP(cp_len+1:end, :);
rx_blocks = fft(rx_ACOOFDMIM);

% MMSE Equalization
rx_blocks = rx_blocks .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);

% ACO-OFDM-IM detection
rx_bits_ACOIM = [];
for sym = 1:numSymbols
    for g = 1:effectiveGroupsACO
        startSubcarrier = 2 + (g-1)*group_size;
        groupSubcarriers = startSubcarrier:startSubcarrier+group_size-1;
        rx_group = rx_blocks(groupSubcarriers, sym);
        
        % Find active subcarriers
        [~, sorted_idx] = sort(abs(rx_group).^2, 'descend');
        detected_active = sort(sorted_idx(1:active_per_group))';
        
        % Find matching combination
        diff = abs(combinations - detected_active);
        [~, comb_idx] = min(sum(diff, 2));
        
        % Decode index bits
        rx_index_bits = de2bi(comb_idx-1, log2(numCombinations), 'left-msb');
        
        % Decode data symbols
        rx_symbols = rx_group(combinations(comb_idx, :));
        rx_data_bits = int2bits(qamdemod(rx_symbols, M, 'UnitAveragePower', true), k);
        
        rx_bits_ACOIM = [rx_bits_ACOIM; rx_index_bits(:); rx_data_bits(:)];
    end
end

% Remove padding if any
if padBits > 0
    rx_bits_ACOIM = rx_bits_ACOIM(1:end-(bitsPerGroupACO-padBits));
end

ber_ACOOFDMIM(snrIdx) = mean(bits_ACOOFDMIM ~= rx_bits_ACOIM(1:length(bits_ACOOFDMIM)));
    %% =================== SCMA ===================
    bits_SCMA = randi([0 1], scma_users*numSymbols*k, 1);
    bits_reshaped = reshape(bits_SCMA, k*scma_users, numSymbols).';
    cp_len_scma = maxDelay;  % Cyclic prefix length
    % SCMA encoding
    tx_SCMA = zeros(scma_resources, numSymbols);
    for sym = 1:numSymbols
        tx_SCMA(1, sym) = scma_codebooks{1}(bi2de(bits_reshaped(sym, 1:k), 'left-msb')+1);
        tx_SCMA(3, sym) = scma_codebooks{2}(bi2de(bits_reshaped(sym, k+1:end), 'left-msb')+1);
    end
    
    % Channel and noise
    rx_SCMA = zeros(scma_resources, numSymbols);
    for sym = 1:numSymbols
        rx_SCMA(:, sym) = filter(channelTaps, 1, tx_SCMA(:, sym));
    end
    rx_SCMA = awgn(rx_SCMA, snr, 'measured');
    
    % SCMA decoding (ML)
    rx_bits_SCMA = zeros(size(bits_SCMA));
    for sym = 1:numSymbols
        min_dist = inf;
        best_bits = zeros(1, k*scma_users);
        
        for b1 = 0:M-1
            for b2 = 0:M-1
                % Encode both users
                cw = zeros(scma_resources, 1);
                cw(1) = scma_codebooks{1}(b1+1);
                cw(3) = scma_codebooks{2}(b2+1);
                
                % Apply channel
                cw_channel = filter(channelTaps, 1, cw);
                
                % Calculate distance
                dist = sum(abs(rx_SCMA(:, sym) - cw_channel).^2);
                
                if dist < min_dist
                    min_dist = dist;
                    best_bits = [de2bi(b1, k, 'left-msb'), de2bi(b2, k, 'left-msb')];
                end
            end
        end
        rx_bits_SCMA((sym-1)*k*scma_users+1:sym*k*scma_users) = best_bits;
    end
    ber_SCMA(snrIdx) = mean(bits_SCMA ~= rx_bits_SCMA);


%% =================== Optical SCMA (DCO-SCMA and ACO-SCMA) ===================
bits_SCMA = randi([0 1], scma_users*numSymbols*k, 1);
bits_reshaped = reshape(bits_SCMA, k*scma_users, numSymbols).';

% Using real-valued constellations with Hermitian symmetry
scma_codebooks_optical = {
    [1, -1, 0.5, -0.5];  % User 1 (real-valued symbols)
    [0.5, -0.5, 1, -1];   % User 2 (real-valued symbols)
};

% SCMA encoding with Hermitian symmetry
   %% DCO-SCMA (DC-Biased Optical SCMA)
    tx_DCOSCMA = zeros(N, numSymbols);
    for sym = 1:numSymbols
        codeword = zeros(N, 1);
        
        % User 1 (resources 1 & 3)
        user1_sym = scma_codebooks_optical{1}(bi2de(bits_reshaped(sym, 1:k), 'left-msb') + 1);
        codeword(1) = user1_sym;
        codeword(3) = user1_sym;
        
        % User 2 (resources 2 & 4)
        user2_sym = scma_codebooks_optical{2}(bi2de(bits_reshaped(sym, k+1:end), 'left-msb') + 1);
        codeword(2) = user2_sym;
        codeword(4) = user2_sym;
        
        % Hermitian symmetry for real signal
        codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
        
        % IFFT and DC bias
        tx_DCOSCMA(:, sym) = ifft(codeword, N, 'symmetric');
    end
    
    % Add DC bias (make signal non-negative)
    dcBias = -min(tx_DCOSCMA(:));
    tx_DCOSCMA = tx_DCOSCMA + dcBias;
    
    % Add CP and transmit
    tx_DCOSCMA_CP = [tx_DCOSCMA(end-cp_len_scma+1:end, :); tx_DCOSCMA];
    rx_DCOSCMA_CP = zeros(size(tx_DCOSCMA_CP));
    for sym = 1:numSymbols
        rx_DCOSCMA_CP(:, sym) = filter(channelTaps, 1, tx_DCOSCMA_CP(:, sym));
    end
    rx_DCOSCMA_CP = awgn(rx_DCOSCMA_CP, snr, 'measured');
    
    % Remove CP and FFT
    rx_DCOSCMA = rx_DCOSCMA_CP(cp_len_scma+1:end, :);
    rx_DCOSCMA_freq = fft(rx_DCOSCMA);
    
    % MMSE Equalization
    H_est = fft(channelTaps, N).';
    rx_DCOSCMA_freq = rx_DCOSCMA_freq .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    
    % ML Detection (DCO-SCMA)
    rx_bits_DCOSCMA = zeros(size(bits_SCMA));
    for sym = 1:numSymbols
        min_dist = inf;
        best_bits = zeros(1, k * scma_users);
        
        for b1 = 0:M-1
            for b2 = 0:M-1
                % Generate codeword
                codeword = zeros(N, 1);
                codeword(1) = scma_codebooks_optical{1}(b1 + 1);
                codeword(3) = scma_codebooks_optical{1}(b1 + 1);
                codeword(2) = scma_codebooks_optical{2}(b2 + 1);
                codeword(4) = scma_codebooks_optical{2}(b2 + 1);
                
                % Hermitian symmetry
                codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
                
                % Modulate and add DC bias
                tx_sym = ifft(codeword, N, 'symmetric') + dcBias;
                
                % Apply channel
                rx_sym = filter(channelTaps, 1, tx_sym);
                
                % Compare with received signal
                dist = sum(abs(rx_DCOSCMA(:, sym) - rx_sym).^2);
                
                if dist < min_dist
                    min_dist = dist;
                    best_bits = [de2bi(b1, k, 'left-msb'), de2bi(b2, k, 'left-msb')];
                end
            end
        end
        rx_bits_DCOSCMA((sym-1)*k*scma_users+1 : sym*k*scma_users) = best_bits;
    end
    ber_DCOSCMA(snrIdx) = mean(bits_SCMA ~= rx_bits_DCOSCMA);
    
    %% ACO-SCMA 
    tx_ACOSCMA = zeros(N, numSymbols);
    for sym = 1:numSymbols
        codeword = zeros(N, 1);
        
        % User 1 (resources 1 & 3)
        user1_sym = scma_codebooks_optical{1}(bi2de(bits_reshaped(sym, 1:k), 'left-msb') + 1);
        codeword(1) = user1_sym;
        codeword(3) = user1_sym;
        
        % User 2 (resources 2 & 4)
        user2_sym = scma_codebooks_optical{2}(bi2de(bits_reshaped(sym, k+1:end), 'left-msb') + 1);
        codeword(2) = user2_sym;
        codeword(4) = user2_sym;
        
        % Hermitian symmetry for real signal
        codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
        
        % IFFT and clip negatives
        tx_ACOSCMA(:, sym) = ifft(codeword, N, 'symmetric');
        tx_ACOSCMA(tx_ACOSCMA < 0) = 0; % ACO clipping
    end
    
    % Add CP and transmit
    tx_ACOSCMA_CP = [tx_ACOSCMA(end-cp_len_scma+1:end, :); tx_ACOSCMA]
    rx_ACOSCMA_CP = zeros(size(tx_ACOSCMA_CP));
    for sym = 1:numSymbols
        rx_ACOSCMA_CP(:, sym) = filter(channelTaps, 1, tx_ACOSCMA_CP(:, sym));
    end
    rx_ACOSCMA_CP = awgn(rx_ACOSCMA_CP, snr, 'measured');
    
    % Remove CP and FFT
    rx_ACOSCMA = rx_ACOSCMA_CP(cp_len_scma+1:end, :);
    rx_ACOSCMA_freq = fft(rx_ACOSCMA);
    
    % MMSE Equalization
    rx_ACOSCMA_freq = rx_ACOSCMA_freq .* conj(H_est) ./ (abs(H_est).^2 + noiseVar);
    
    % ML Detection (ACO-SCMA)
    rx_bits_ACOSCMA = zeros(size(bits_SCMA));
    for sym = 1:numSymbols
        min_dist = inf;
        best_bits = zeros(1, k * scma_users);
        
        for b1 = 0:M-1
            for b2 = 0:M-1
                % Generate codeword
                codeword = zeros(N, 1);
                codeword(1) = scma_codebooks_optical{1}(b1 + 1);
                codeword(3) = scma_codebooks_optical{1}(b1 + 1);
                codeword(2) = scma_codebooks_optical{2}(b2 + 1);
                codeword(4) = scma_codebooks_optical{2}(b2 + 1);
                
                % Hermitian symmetry
                codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
                
                % Modulate and clip negatives
                tx_sym = ifft(codeword, N, 'symmetric');
                tx_sym(tx_sym < 0) = 0;
                
                % Apply channel
                rx_sym = filter(channelTaps, 1, tx_sym);
                
                % Compare with received signal
                dist = sum(abs(rx_ACOSCMA(:, sym) - rx_sym).^2);
                
                if dist < min_dist
                    min_dist = dist;
                    best_bits = [de2bi(b1, k, 'left-msb'), de2bi(b2, k, 'left-msb')];
                end
            end
        end
        rx_bits_ACOSCMA((sym-1)*k*scma_users+1 : sym*k*scma_users) = best_bits;
    end
    ber_ACOSCMA(snrIdx) = mean(bits_SCMA ~= rx_bits_ACOSCMA);
end


%% Plot Results - RF vs Optical Comparison
figure;
% RF waveforms (solid blue lines)
semilogy(SNR_dB, ber_OFDM, '-o', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OFDM'); hold on;
semilogy(SNR_dB, ber_OTFS, '-x', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OTFS');
semilogy(SNR_dB, ber_OCDM, '-+', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OCDM');
semilogy(SNR_dB, ber_OFDMIM, '-s', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OFDM-IM');
semilogy(SNR_dB, ber_SCMA, '-d', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'SCMA');

% Optical waveforms (dashed green lines)
semilogy(SNR_dB, ber_ACOOFDM, '--o', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(SNR_dB, ber_ACOTFS, '--x', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(SNR_dB, ber_ACOCDM, '--+', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(SNR_dB, ber_DCOOFDMIM, '-s', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(SNR_dB, ber_ACOSCMA, '--d', 'Color',[0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('SNR (dB)'); ylabel('BER');
title('RF vs Optical Waveform Comparison');
legend('Location', 'southwest', 'NumColumns', 2);
grid on; ylim([1e-4, 1]);
set(gca, 'FontSize', 12);

%% Plot Results - ACO vs DCO Variants of Optical Waveforms
figure;
% DCO variants (solid green lines)
semilogy(SNR_dB, ber_DCOOFDM, '-o', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM'); hold on;
semilogy(SNR_dB, ber_DCOTFS, '-x', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
semilogy(SNR_dB, ber_DCOCDM, '-+', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(SNR_dB, ber_DCOOFDMIM, '-s', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(SNR_dB, ber_DCOSCMA, '-d', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

% ACO variants (dashed green lines - slightly darker)
semilogy(SNR_dB, ber_ACOOFDM, '--o', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(SNR_dB, ber_ACOTFS, '--x', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(SNR_dB, ber_ACOCDM, '--+', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(SNR_dB, ber_ACOOFDMIM, '--s', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(SNR_dB, ber_ACOSCMA, '--d', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('SNR (dB)'); ylabel('BER');
title('ACO vs DCO Variants of Optical Waveforms');
legend('Location', 'southwest', 'NumColumns', 2);
grid on; ylim([1e-4, 1]);
set(gca, 'FontSize', 12);


%% 2nd version:

%% Plot Results - RF vs Optical Comparison
figure;
% RF waveforms (solid blue lines)
semilogy(SNR_dB, ber_OFDM, '-o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'OFDM'); hold on;
semilogy(SNR_dB, ber_OTFS, '-x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'OTFS');
semilogy(SNR_dB, ber_OCDM, '-+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'OCDM');
semilogy(SNR_dB, ber_OFDMIM, '-s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'OFDM-IM');
semilogy(SNR_dB, ber_SCMA, '-d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'SCMA');

% Optical waveforms (dashed green lines)
semilogy(SNR_dB, ber_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(SNR_dB, ber_ACOTFS, '--x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(SNR_dB, ber_ACOCDM, '--+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(SNR_dB, ber_DCOOFDMIM, '-.s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(SNR_dB, ber_ACOSCMA, '--d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('SNR (dB)'); ylabel('BER');
% title('RF vs Optical Waveform Comparison');
lgd=legend('Location', 'southwest', 'NumColumns', 2);
% lgd.Color = 'none'; % Makes the background itself transparent
set(lgd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.8]));

grid on; ylim([1e-4, 1]);
set(gca, 'FontSize', 12);
ylim([1e-3, 1])

%% Plot Results - ACO vs DCO Variants of Optical Waveforms
figure;
% DCO variants (solid green lines)
semilogy(SNR_dB, ber_DCOOFDM, '-.o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM'); hold on;
semilogy(SNR_dB, ber_DCOTFS, '-.x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
semilogy(SNR_dB, ber_DCOCDM, '-.+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(SNR_dB, ber_DCOOFDMIM, '-.s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(SNR_dB, ber_DCOSCMA, '-.d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

% ACO variants (dashed green lines - slightly darker)
semilogy(SNR_dB, ber_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(SNR_dB, ber_ACOTFS, '--x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(SNR_dB, ber_ACOCDM, '--+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(SNR_dB, ber_ACOOFDMIM, '--s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(SNR_dB, ber_ACOSCMA, '--d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('SNR (dB)'); ylabel('BER');
% title('ACO vs DCO Variants of Optical Waveforms');
lgd=legend('Location', 'southwest', 'NumColumns', 2);
set(lgd.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;0.8]));
grid on; ylim([1e-4, 1]);
set(gca, 'FontSize', 12);
ylim([1e-3, 1])


%% Helper Functions
function int = bits2int(bits, k)
    if size(bits,2) ~= k
        bits = reshape(bits, k, []).';
    end
    int = bi2de(bits, 'left-msb');
end

function bits = int2bits(int, k)
    bits = de2bi(int, k, 'left-msb').';
    bits = bits(:);
end