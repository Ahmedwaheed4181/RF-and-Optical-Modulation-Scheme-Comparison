clc; clear; close all;
rng(0); 
%% Common Parameters
N = 64;                 % Number of subcarriers/resources (power of 2 for Hermitian symmetry)
M = 4;                  % QPSK modulation
k = log2(M);            % Bits per symbol
numSymbols = 100;        % Number of symbols per frame
cp_len = N/8;           % Cyclic prefix length
numPAPRSamples = 1e4; % Number of PAPR samples
ccdfPoints = 20;        % Number of points for CCDF curve

%% OFDM-IM Parameters (for IM variants)
group_size = 8;         % Subcarriers per group
active_per_group = 4;   % Active subcarriers per group
numGroups = N / group_size;
combinations = nchoosek(1:group_size, active_per_group);
numCombinations = size(combinations, 1);
bits_per_group = log2(numCombinations) + active_per_group*k;

%% SCMA Parameters 
scma_users = 2;         % Number of users
scma_resources = 4;     % Number of resources
scma_codebooks = {      % SCMA codebooks
    [1+1i, -1+1i, 1-1i, -1-1i]/sqrt(2);  % User 1
    [1+1i, 1-1i, -1+1i, -1-1i]/sqrt(2)   % User 2
};
scma_codebooks_optical = { % SCMA codebooks for optical (real-valued)
    [1, -1, 0.5, -0.5];    % User 1
    [0.5, -0.5, 1, -1]     % User 2
};

%% Preallocate PAPR arrays
% Conventional waveforms
papr_OFDM = zeros(numPAPRSamples, 1);
papr_OTFS = zeros(numPAPRSamples, 1);
papr_OCDM = zeros(numPAPRSamples, 1);
papr_OFDMIM = zeros(numPAPRSamples, 1);
papr_SCMA = zeros(numPAPRSamples, 1);

% Optical waveforms
papr_DCOOFDM = zeros(numPAPRSamples, 1);
papr_ACOOFDM = zeros(numPAPRSamples, 1);
papr_DCOTFS = zeros(numPAPRSamples, 1);
papr_ACOTFS = zeros(numPAPRSamples, 1);
papr_DCOCDM = zeros(numPAPRSamples, 1);
papr_ACOCDM = zeros(numPAPRSamples, 1);
papr_DCOOFDMIM = zeros(numPAPRSamples, 1);
papr_ACOOFDMIM = zeros(numPAPRSamples, 1);
papr_DCOSCMA = zeros(numPAPRSamples, 1);
papr_ACOSCMA = zeros(numPAPRSamples, 1);

%% Generate PAPR samples
for i = 1:numPAPRSamples
    %% ============= Conventional Waveforms =============
    %% OFDM PAPR
    symbols_OFDM = qammod(randi([0 M-1], N, 1), M, 'UnitAveragePower', true);
    tx_OFDM = ifft(symbols_OFDM);
    papr_OFDM(i) = 10*log10(max(abs(tx_OFDM).^2) / mean(abs(tx_OFDM).^2));
    
    %% OTFS PAPR
    symbols_OTFS = qammod(randi([0 M-1], N, 1), M, 'UnitAveragePower', true);
    dd_grid = reshape(symbols_OTFS, N, 1);
    tx_OTFS = ifft(fft(dd_grid, [], 2), [], 1);
    papr_OTFS(i) = 10*log10(max(abs(tx_OTFS).^2) / mean(abs(tx_OTFS).^2));
    
    %% OCDM PAPR
    symbols_OCDM = qammod(randi([0 M-1], N, 1), M, 'UnitAveragePower', true);
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N);
    tx_OCDM = ifft(symbols_OCDM .* chirp);
    papr_OCDM(i) = 10*log10(max(abs(tx_OCDM).^2) / mean(abs(tx_OCDM).^2));
    
    %% OFDM-IM PAPR
    index_decimal = randi([1 numCombinations], numGroups, 1);
    mod_data = qammod(randi([0 M-1], numGroups*active_per_group, 1), M, 'UnitAveragePower', true);
    mod_data = reshape(mod_data, numGroups, active_per_group);
    
    tx_blocks = zeros(N, 1);
    for g = 1:numGroups
        idx = (g-1)*group_size + combinations(index_decimal(g), :);
        tx_blocks(idx) = mod_data(g, :);
    end
    tx_OFDMIM = ifft(tx_blocks, N);
    papr_OFDMIM(i) = 10*log10(max(abs(tx_OFDMIM).^2) / mean(abs(tx_OFDMIM).^2));
    
    %% SCMA PAPR
    tx_SCMA_signal = [];
    scma_symbols_per_papr = 100; % SCMA symbols per PAPR calculation
  
cp_len = min(8, scma_resources-1); % Ensure CP length is valid for SCMA
    for sym = 1:scma_symbols_per_papr
        % Create random SCMA codeword
        cw = zeros(scma_resources, 1);
        active_resources = randperm(scma_resources, 2); % 2 active resources
        cw(active_resources(1)) = scma_codebooks{1}(randi([1 M]));
        cw(active_resources(2)) = scma_codebooks{2}(randi([1 M]));
        
        % Convert to time domain with CP
        tx_sym = ifft(cw);
        % Ensure CP length doesn't exceed symbol length
        cp_len_actual = min(cp_len, length(tx_sym)-1);
        tx_sym_cp = [tx_sym(end-cp_len_actual+1:end); tx_sym];
        tx_SCMA_signal = [tx_SCMA_signal; tx_sym_cp];
    end
    papr_SCMA(i) = 10*log10(max(abs(tx_SCMA_signal).^2) / mean(abs(tx_SCMA_signal).^2));
    
    %% ============= Optical Waveforms =============
    %% =================== DCO-OFDM ===================
    symbols_DCOOFDM = qammod(randi([0 M-1], (N/2-1)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_DCOOFDM, (N/2-1), numSymbols);
    data2 = [zeros(1,numSymbols); data; zeros(1,numSymbols); flipud(conj(data))];
    tx_DCOOFDM = ifft(data2);
    dcBias = -min(real(tx_DCOOFDM(:)));
    tx_DCOOFDM = real(tx_DCOOFDM) + dcBias;
    papr_DCOOFDM(i) = 10*log10(max(abs(tx_DCOOFDM(:)).^2) / mean(abs(tx_DCOOFDM(:)).^2));
    
    %% =================== ACO-OFDM ===================
    symbols_ACOOFDM = qammod(randi([0 M-1], (N/4)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_ACOOFDM, (N/4), numSymbols);
    data2 = zeros(N,numSymbols);
    for kk=1:N/2
        if ~rem(kk,2)
            data2(kk,:) = data(kk/2,:);
        end
    end
    data2(N/2+2:end,:) = flipud(conj(data2(2:N/2,:)));
    tx_ACOOFDM = ifft(data2);
    tx_ACOOFDM(tx_ACOOFDM<0) = 0;
    papr_ACOOFDM(i) = 10*log10(max(abs(tx_ACOOFDM(:)).^2) / mean(abs(tx_ACOOFDM(:)).^2));
    
    %% =================== DCO-OTFS ===================
    symbols_DCOTFS = qammod(randi([0 M-1], (N/2-1)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_DCOTFS, (N/2-1), numSymbols);
    data2 = [zeros(1,numSymbols); data; zeros(1,numSymbols); flipud(conj(data))];
    dd_grid = zeros(N, numSymbols);
    dd_grid(2:N/2,:) = data2(2:N/2,:);
    dd_grid(N/2+2:end,:) = conj(flipud(data2(2:N/2,:)));
    tx_DCOTFS = ifft(fft(dd_grid, [], 2), [], 1, 'symmetric');
    dcBias = -min(tx_DCOTFS(:));
    tx_DCOTFS = tx_DCOTFS + dcBias;
    papr_DCOTFS(i) = 10*log10(max(abs(tx_DCOTFS(:)).^2) / mean(abs(tx_DCOTFS(:)).^2));
    
    %% =================== ACO-OTFS ===================
    symbols_ACOTFS = qammod(randi([0 M-1], (N/4)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_ACOTFS, (N/4), numSymbols);
    dd_grid = zeros(N, numSymbols);
    for kk = 1:N/2
        if mod(kk,2) == 0
            dd_grid(kk,:) = data(kk/2,:);
        end
    end
    dd_grid(N/2+2:end,:) = conj(flipud(dd_grid(2:N/2,:)));
    tx_ACOTFS = ifft(fft(dd_grid, [], 2), [], 1, 'symmetric');
    tx_ACOTFS(tx_ACOTFS < 0) = 0;
    papr_ACOTFS(i) = 10*log10(max(abs(tx_ACOTFS(:)).^2) / mean(abs(tx_ACOTFS(:)).^2));
    
    %% =================== DCO-OCDM ===================
    symbols_DCOCDM = qammod(randi([0 M-1], (N/2-1)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_DCOCDM, (N/2-1), numSymbols);
    data2 = [zeros(1,numSymbols); data; zeros(1,numSymbols); flipud(conj(data))];
    input_sym = zeros(N, numSymbols);
    input_sym(2:N/2,:) = data2(2:N/2,:);
    input_sym(N/2+2:end,:) = conj(flipud(data2(2:N/2,:)));
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N);
    tx_DCOCDM = ifft(input_sym .* chirp, [], 1, 'symmetric');
    dcBias = -min(tx_DCOCDM(:));
    tx_DCOCDM = tx_DCOCDM + dcBias;
    papr_DCOCDM(i) = 10*log10(max(abs(tx_DCOCDM(:)).^2) / mean(abs(tx_DCOCDM(:)).^2));
    
    %% =================== ACO-OCDM ===================
    symbols_ACOCDM = qammod(randi([0 M-1], (N/4)*numSymbols, 1), M, 'UnitAveragePower', true);
    data = reshape(symbols_ACOCDM, (N/4), numSymbols);
    input_sym = zeros(N, numSymbols);
    for kk = 1:N/2
        if mod(kk,2) == 0
            input_sym(kk,:) = data(kk/2,:);
        end
    end
    input_sym(N/2+2:end,:) = conj(flipud(input_sym(2:N/2,:)));
    n = (0:N-1).';
    chirp = exp(1j*pi*(n.^2)/N);
    tx_ACOCDM = ifft(input_sym .* chirp, [], 1, 'symmetric');
    tx_ACOCDM(tx_ACOCDM < 0) = 0;
    papr_ACOCDM(i) = 10*log10(max(abs(tx_ACOCDM(:)).^2) / mean(abs(tx_ACOCDM(:)).^2));
    
    %% =================== DCO-OFDM-IM ===================
    effective_N = N/2-1;
    effective_numGroups = floor(effective_N/group_size);
    index_decimal = randi([1 numCombinations], effective_numGroups*numSymbols, 1);
    mod_data = qammod(randi([0 M-1], effective_numGroups*numSymbols*active_per_group, 1), M, 'UnitAveragePower', true);
    mod_data = reshape(mod_data, effective_numGroups*numSymbols, active_per_group);
    
    tx_blocks = zeros(N, numSymbols);
    for sym = 1:numSymbols
        for g = 1:effective_numGroups
            idx = (g-1)*group_size + combinations(index_decimal((sym-1)*effective_numGroups + g), :) + 1;
            tx_blocks(idx, sym) = mod_data((sym-1)*effective_numGroups + g, :);
        end
        tx_blocks(N/2+2:end, sym) = conj(flipud(tx_blocks(2:N/2, sym)));
    end
    tx_DCOOFDMIM = ifft(tx_blocks, N, 'symmetric');
    dcBias = -min(tx_DCOOFDMIM(:));
    tx_DCOOFDMIM = tx_DCOOFDMIM + dcBias;
    papr_DCOOFDMIM(i) = 10*log10(max(abs(tx_DCOOFDMIM(:)).^2) / mean(abs(tx_DCOOFDMIM(:)).^2));
    
    %% =================== ACO-OFDM-IM ===================
    numEvenSubcarriers = N/2;
    effectiveGroupsACO = floor(numEvenSubcarriers/group_size);
    index_decimal = randi([1 numCombinations], effectiveGroupsACO*numSymbols, 1);
    mod_data = qammod(randi([0 M-1], effectiveGroupsACO*numSymbols*active_per_group, 1), M, 'UnitAveragePower', true);
    mod_data = reshape(mod_data, effectiveGroupsACO*numSymbols, active_per_group);
    
    tx_blocks = zeros(N, numSymbols);
    for sym = 1:numSymbols
        for g = 1:effectiveGroupsACO
            startSubcarrier = 2 + (g-1)*group_size;
            activeSubcarriers = startSubcarrier + combinations(index_decimal((sym-1)*effectiveGroupsACO + g), :) - 1;
            tx_blocks(activeSubcarriers, sym) = mod_data((sym-1)*effectiveGroupsACO + g, :);
        end
        tx_blocks(N/2+2:end, sym) = conj(flipud(tx_blocks(2:N/2, sym)));
    end
    tx_ACOOFDMIM = ifft(tx_blocks, N, 'symmetric');
    tx_ACOOFDMIM(tx_ACOOFDMIM < 0) = 0;
    papr_ACOOFDMIM(i) = 10*log10(max(abs(tx_ACOOFDMIM(:)).^2) / mean(abs(tx_ACOOFDMIM(:)).^2));
    
    %% =================== DCO-SCMA ===================
    tx_DCOSCMA = zeros(N, numSymbols);
    for sym = 1:numSymbols
        codeword = zeros(N, 1);
        user1_sym = scma_codebooks_optical{1}(randi([1 M]));
        codeword(1) = user1_sym;
        codeword(3) = user1_sym;
        user2_sym = scma_codebooks_optical{2}(randi([1 M]));
        codeword(2) = user2_sym;
        codeword(4) = user2_sym;
        codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
        tx_DCOSCMA(:, sym) = ifft(codeword, N, 'symmetric');
    end
    dcBias = -min(tx_DCOSCMA(:));
    tx_DCOSCMA = tx_DCOSCMA + dcBias;
    papr_DCOSCMA(i) = 10*log10(max(abs(tx_DCOSCMA(:)).^2) / mean(abs(tx_DCOSCMA(:)).^2));
    
    %% =================== ACO-SCMA ===================
    tx_ACOSCMA = zeros(N, numSymbols);
    for sym = 1:numSymbols
        codeword = zeros(N, 1);
        user1_sym = scma_codebooks_optical{1}(randi([1 M]));
        codeword(1) = user1_sym;
        codeword(3) = user1_sym;
        user2_sym = scma_codebooks_optical{2}(randi([1 M]));
        codeword(2) = user2_sym;
        codeword(4) = user2_sym;
        codeword(N/2+2:end) = conj(flipud(codeword(2:N/2)));
        tx_ACOSCMA(:, sym) = ifft(codeword, N, 'symmetric');
    end
    tx_ACOSCMA(tx_ACOSCMA < 0) = 0;
    papr_ACOSCMA(i) = 10*log10(max(abs(tx_ACOSCMA(:)).^2) / mean(abs(tx_ACOSCMA(:)).^2));
end

%% Calculate CCDF
[paprAxis, ccdf_OFDM] = calculateCCDF(papr_OFDM, ccdfPoints);
[~, ccdf_OTFS] = calculateCCDF(papr_OTFS, ccdfPoints);
[~, ccdf_OCDM] = calculateCCDF(papr_OCDM, ccdfPoints);
[~, ccdf_OFDMIM] = calculateCCDF(papr_OFDMIM, ccdfPoints);
[~, ccdf_SCMA] = calculateCCDF(papr_SCMA, ccdfPoints);

[~, ccdf_DCOOFDM] = calculateCCDF(papr_DCOOFDM, ccdfPoints);
[~, ccdf_ACOOFDM] = calculateCCDF(papr_ACOOFDM, ccdfPoints);
[~, ccdf_DCOTFS] = calculateCCDF(papr_DCOTFS, ccdfPoints);
[~, ccdf_ACOTFS] = calculateCCDF(papr_ACOTFS, ccdfPoints);
[~, ccdf_DCOCDM] = calculateCCDF(papr_DCOCDM, ccdfPoints);
[~, ccdf_ACOCDM] = calculateCCDF(papr_ACOCDM, ccdfPoints);
[~, ccdf_DCOOFDMIM] = calculateCCDF(papr_DCOOFDMIM, ccdfPoints);
[~, ccdf_ACOOFDMIM] = calculateCCDF(papr_ACOOFDMIM, ccdfPoints);
[~, ccdf_DCOSCMA] = calculateCCDF(papr_DCOSCMA, ccdfPoints);
[~, ccdf_ACOSCMA] = calculateCCDF(papr_ACOSCMA, ccdfPoints);

%% Plot CCDF for Different Waveforms (using default MATLAB colors)
figure;

% RF waveforms (solid blue lines with different markers)
semilogy(paprAxis, ccdf_OFDM, '-o', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OFDM'); hold on;
semilogy(paprAxis, ccdf_ACOOFDM, '--o', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(paprAxis, ccdf_OTFS, '-x', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OTFS');
semilogy(paprAxis, ccdf_OCDM, '-+', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OCDM');
semilogy(paprAxis, ccdf_OFDMIM, '-s', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'OFDM-IM');
semilogy(paprAxis, ccdf_SCMA, '-d', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'SCMA');
semilogy(paprAxis, ccdf_ACOTFS, '--x', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(paprAxis, ccdf_DCOCDM, '-+', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(paprAxis, ccdf_ACOOFDMIM, '--s', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(paprAxis, ccdf_DCOSCMA, '-d', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

xlabel('PAPR (dB)'); ylabel('CCDF (Pr(PAPR > x))');
title('PAPR Comparison of RF Waveforms');
legend('Location', 'northeast');
grid on;
ylim([1e-3 1]);
xlim([4 12]);

%% Plot CCDF for Optical Waveforms
figure;

% DCO variants (solid green lines with different markers)
semilogy(paprAxis, ccdf_DCOOFDM, '-o', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM'); hold on;
semilogy(paprAxis, ccdf_DCOTFS, '-x', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
semilogy(paprAxis, ccdf_DCOCDM, '-+', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(paprAxis, ccdf_DCOOFDMIM, '-s', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(paprAxis, ccdf_DCOSCMA, '-d', 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

% ACO variants (dashed green lines with different markers)
semilogy(paprAxis, ccdf_ACOOFDM, '--o', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(paprAxis, ccdf_ACOTFS, '--x', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(paprAxis, ccdf_ACOCDM, '--+', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(paprAxis, ccdf_ACOOFDMIM, '--s', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(paprAxis, ccdf_ACOSCMA, '--d', 'Color', [0.1960, 0.4740, 0.0880], 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('PAPR (dB)'); ylabel('CCDF (Pr(PAPR > x))');
title('PAPR Comparison of Optical Waveforms');
legend('Location', 'northeast');
grid on;
ylim([1e-4 1]);
xlim([4 16]);


%% 2nd Version:
%% Plot CCDF for Different Waveforms (using default MATLAB colors)
figure;

% RF waveforms (solid blue lines with different markers)
semilogy(paprAxis, ccdf_OFDM, '-o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'OFDM'); hold on;
semilogy(paprAxis, ccdf_OTFS, '-x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'OTFS');
semilogy(paprAxis, ccdf_OCDM, '-+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'OCDM');
semilogy(paprAxis, ccdf_OFDMIM, '-s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'OFDM-IM');
semilogy(paprAxis, ccdf_SCMA, '-d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'SCMA');

semilogy(paprAxis, ccdf_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(paprAxis, ccdf_DCOTFS, '-.x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
semilogy(paprAxis, ccdf_DCOCDM, '-.+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(paprAxis, ccdf_ACOOFDMIM, '-.s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(paprAxis, ccdf_DCOSCMA, '-.d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

xlabel('PAPR (dB)'); ylabel('CCDF');
% title('PAPR Comparison of RF Waveforms');
legend('Location', 'northeast', 'NumColumns', 2);
grid on;
ylim([1e-3 1]);
xlim([6 11]);

%% Plot CCDF for Optical Waveforms
figure;

% DCO variants (solid green lines with different markers)
semilogy(paprAxis, ccdf_DCOOFDM, '-.o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM'); hold on;
semilogy(paprAxis, ccdf_DCOTFS, '-.x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'DCO-OTFS');
semilogy(paprAxis, ccdf_DCOCDM, '-.+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'DCO-OCDM');
semilogy(paprAxis, ccdf_DCOOFDMIM, '-.s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'DCO-OFDM-IM');
semilogy(paprAxis, ccdf_DCOSCMA, '-.d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'DCO-SCMA');

% ACO variants (dashed green lines with different markers)
semilogy(paprAxis, ccdf_ACOOFDM, '--o', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM');
semilogy(paprAxis, ccdf_ACOTFS, '--x', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'ACO-OTFS');
semilogy(paprAxis, ccdf_ACOCDM, '--+', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ACO-OCDM');
semilogy(paprAxis, ccdf_ACOOFDMIM, '--s', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'ACO-OFDM-IM');
semilogy(paprAxis, ccdf_ACOSCMA, '--d', 'Color', 'c', 'LineWidth', 1.5, 'DisplayName', 'ACO-SCMA');

xlabel('PAPR (dB)'); ylabel('CCDF');
% title('PAPR Comparison of Optical Waveforms');
legend('Location', 'northeast', 'NumColumns', 2);
grid on;
ylim([1e-3 1]);
xlim([6 11]);


%% Helper Functions
function [paprAxis, ccdf] = calculateCCDF(paprSamples, numPoints)
    [counts, edges] = histcounts(paprSamples, numPoints, 'Normalization', 'cdf');
    paprAxis = edges(1:end-1) + diff(edges)/2;
    ccdf = 1 - counts;
end