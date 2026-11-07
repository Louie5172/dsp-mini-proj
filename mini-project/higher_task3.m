%% Task 3 - Remove Fundamental Frequency

% Original speech
[s_original, fs_recorded] = audioread('Female16K.wav');

assert(fs_recorded == fs, ...
    'Speech sampling frequency must match the cochlear model.');

% Design elliptic bandpass filter: 320 Hz to 7500 Hz
Fstop = [280 7800];
Fpass = [320 7500];

Rp = 0.1;  % Passband ripple (dB)
Rs = 80;   % Stopband attenuation (dB)

[filter_order, Wn] = ellipord( ...
    Fpass/(fs/2), Fstop/(fs/2), Rp, Rs);

[z_bp, p_bp, k_bp] = ellip( ...
    filter_order, Rp, Rs, Wn, 'bandpass');

[sos_bp, g_bp] = zp2sos(z_bp, p_bp, k_bp);

fprintf('Task 3 elliptic prototype order: %d\n', filter_order);

% Apply filter
s_filtered = sosfilt(sos_bp, g_bp*s_original);


%% Task 3 - Plot Bandpass Filter Response

[H_bp, f_bp] = freqz(sos_bp, 16384, fs);
H_bp = g_bp * H_bp;

figure;

plot(f_bp, 20*log10(max(abs(H_bp), eps)));
grid on;

xline(280, '--r', '280 Hz Stopband');
xline(320, '--g', '320 Hz Passband');
xline(7500, '--g', '7500 Hz Passband');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Elliptic Bandpass Filter Response');

xlim([0 8000]);
ylim([-120 5]);


%% Task 3 - Listen to Original and Filtered Speech

disp('Playing original speech...');
soundsc(s_original, fs);
pause(length(s_original)/fs + 1);

disp('Playing filtered speech...');
soundsc(s_filtered, fs);


%% Task 3 - Compare Frequency Spectra

N_fft = length(s_original);

Y_original = abs(fft(s_original));
Y_filtered = abs(fft(s_filtered));

f = (0:floor(N_fft/2)) * fs/N_fft;

figure;

subplot(2,1,1);
plot(f, Y_original(1:length(f)));
grid on;
xlim([0 1500]);
title('Original Speech Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2,1,2);
plot(f, Y_filtered(1:length(f)));
xline(320, '--r', '320 Hz Passband Edge');
grid on;
xlim([0 1500]);
title('Elliptic Bandpass Filtered Speech Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%% Task 3 - Run Filtered Speech Through Cochlear Filter Bank

N = length(s_filtered);

displacement_filtered = zeros(N, n);
pressure_filtered = zeros(N, n);

y = s_filtered;

for i = 1:n

    b = coeff(i, Gz_idx) * ...
        [1, -coeff(i, a1_idx), coeff(i, a2_idx)];

    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement_filtered(:, i) = y;

    y = filter(b, 1, y);
    pressure_filtered(:, i) = y;
end


%% Task 3 - Spatial Differentiation

spatial_filtered1 = displacement_filtered(:, 1:end-1) ...
                  - displacement_filtered(:, 2:end);

spatial_filtered2 = spatial_filtered1(:, 1:end-1) ...
                  - spatial_filtered1(:, 2:end);


%% Task 3 - Inner Hair Cell

fc = 30;
c_0 = exp(-2*pi*fc/fs);

e_filtered = max(spatial_filtered2, 0);

energy_filtered = zeros(size(e_filtered));

energy_filtered(1, :) = (1-c_0)*e_filtered(1, :);

for i = 2:N
    energy_filtered(i, :) = (1-c_0)*e_filtered(i, :) ...
                          + c_0*energy_filtered(i-1, :);
end

% Read at the end of each 20 ms frame
num_frames_filtered = floor(N/frame_size);
frame_ends_filtered = (1:num_frames_filtered)*frame_size;

energy_frames_filtered = energy_filtered(frame_ends_filtered, :);


%% Task 3 - HPS Pitch Estimation

hps_filtered = zeros(num_frames_filtered, length(fundamental_freqs));
pitch_filtered = zeros(num_frames_filtered, 1);

for frame_idx = 1:num_frames_filtered

    energy_frame = energy_frames_filtered(frame_idx, :);

    hps = prod(energy_frame(harmonic_filters), 2);

    hps_filtered(frame_idx, :) = hps';

    [~, peak_idx] = max(hps);

    pitch_filtered(frame_idx) = fundamental_freqs(peak_idx);
end


%% Task 3 - Compare Inner Hair Cell Outputs

frames_to_test = [20, 30];

for j = 1:length(frames_to_test)

    frame_num = frames_to_test(j);

    if frame_num > min(size(energy_frames,1), num_frames_filtered)
        continue;
    end

    % Original pitch identified in Task 2
    [~, row_idx] = min(abs(fundamental_freqs - pitch(frame_num)));
    pitch_filter = fundamental_filters(row_idx);

    figure;

    subplot(3,1,1);
    plot(energy_frames(frame_num, :));
    xline(pitch_filter, '--r', 'Fundamental');
    grid on;
    xlabel('Filter Number');
    ylabel('Hair Cell Output');
    title(sprintf('Original Speech - Frame %d', frame_num));

    subplot(3,1,2);
    plot(energy_frames_filtered(frame_num, :));
    xline(pitch_filter, '--r', 'Fundamental');
    grid on;
    xlabel('Filter Number');
    ylabel('Hair Cell Output');
    title(sprintf('Filtered Speech - Frame %d', frame_num));

    subplot(3,1,3);
    plot(fundamental_filters, hps_filtered(frame_num, :));
    grid on;
    xlabel('Filter Number');
    ylabel('HPS');
    title(sprintf('Filtered Speech HPS - Frame %d', frame_num));

    fprintf('\nFrame %d\n', frame_num);
    fprintf('Original pitch: %.2f Hz\n', pitch(frame_num));
    fprintf('Filtered pitch: %.2f Hz\n', pitch_filtered(frame_num));
    fprintf('Fundamental filter: %d\n', pitch_filter);
    fprintf('Original energy: %.6g\n', ...
        energy_frames(frame_num, pitch_filter));
    fprintf('Filtered energy: %.6g\n', ...
        energy_frames_filtered(frame_num, pitch_filter));

end


%% Task 3 - Compare Pitch Estimates Across Speech

figure;

plot(1:length(pitch), pitch, '-o', ...
    'DisplayName', 'Original Speech');

hold on;

plot(1:length(pitch_filtered), pitch_filtered, '-o', ...
    'DisplayName', 'Filtered Speech');

hold off;
grid on;
legend;

xlabel('Frame Number');
ylabel('Estimated Fundamental Frequency (Hz)');
title('Pitch Estimation With and Without Fundamental');
ylim([80 420]);


%% Task 4 - Regenerate Missing Fundamental

% B: Already filtered speech
B = s_filtered(:);

% C: Nonlinear operation
C = B.^2;

% D: Elliptic band-pass filter (60-400 Hz)

Fstop_regen = [40 450];
Fpass_regen = [60 400];

Rp_regen = 0.1;
Rs_regen = 60;

[order_regen, Wn_regen] = ellipord( ...
    Fpass_regen/(fs/2), ...
    Fstop_regen/(fs/2), ...
    Rp_regen, Rs_regen);

[z_regen, p_regen, k_regen] = ellip( ...
    order_regen, Rp_regen, Rs_regen, ...
    Wn_regen, 'bandpass');

[sos_regen, g_regen] = zp2sos( ...
    z_regen, p_regen, k_regen);

fprintf('Task 4 elliptic prototype order: %d\n', order_regen);

D = sosfilt(sos_regen, g_regen*C);

% definining a scaling factor, so that this is still picked up as best one
% by the HPS
alpha = 7;

% E: Recombine regenerated fundamental with filtered speech
E = B + alpha * D;


%% Task 4 - Plot Regeneration Filter Response

[H_regen, f_regen] = freqz(sos_regen, 16384, fs);
H_regen = g_regen * H_regen;

figure;

plot(f_regen, 20*log10(max(abs(H_regen), eps)));
grid on;

xline(60, '--g', '60 Hz');
xline(400, '--g', '400 Hz');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Elliptic Regeneration Bandpass Filter');

xlim([0 1000]);
ylim([-120 5]);


%% Task 4 - FFT at Points A, B, C, D and E

signals = {s_original(:), B, C, D, E};

names = {'A - Original', ...
         'B - Filtered', ...
         'C - Nonlinear', ...
         'D - Regenerated', ...
         'E - Recombined'};

% Inspect the same voiced frame at every point
frame_num = 20;
idx = (frame_num-1)*frame_size + (1:frame_size);

N_fft = 8192;
f = (0:N_fft/2)*fs/N_fft;

figure;

for i = 1:5
    subplot(5,1,i);

    Y = abs(fft(signals{i}(idx), N_fft));

    plot(f, Y(1:N_fft/2+1));
    xlim([0 1500]);
    grid on;

    title(names{i});
    xlabel('Frequency (Hz)');
    ylabel('FFT Magnitude');
end


%% Task 4 - Run E Through Cochlear Filter Bank

N_E = length(E);

displacement_E = zeros(N_E, n);

y = E;

for i = 1:n

    b = coeff(i, Gz_idx) * ...
        [1, -coeff(i, a1_idx), coeff(i, a2_idx)];

    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement_E(:, i) = y;

    y = filter(b, 1, y);
end


%% Task 4 - Spatial Differentiation and Hair Cell

spatial_E = diff(displacement_E, 2, 2);

% Half-wave rectification and low-pass filtering
e_E = max(spatial_E, 0);
energy_E = filter(1-c_0, [1, -c_0], e_E);

% Read at the end of each 20 ms frame
frame_ends_E = frame_size:frame_size:N_E;
energy_frames_E = energy_E(frame_ends_E, :);


%% Task 4 - HPS Pitch Estimation at E

num_frames_E = length(frame_ends_E);

hps_E = zeros(num_frames_E, length(fundamental_freqs));
pitch_E = zeros(num_frames_E, 1);

for frame_idx = 1:num_frames_E

    energy_frame = energy_frames_E(frame_idx, :);

    hps = prod(energy_frame(harmonic_filters), 2);

    hps_E(frame_idx, :) = hps';

    [~, peak_idx] = max(hps);

    pitch_E(frame_idx) = fundamental_freqs(peak_idx);
end

% Pitch period
pitch_period_E = fs ./ pitch_E;

fprintf('Frame %d:\n', frame_num);
fprintf('Original pitch: %.2f Hz\n', pitch(frame_num));
fprintf('Filtered pitch: %.2f Hz\n', pitch_filtered(frame_num));
fprintf('Regenerated pitch: %.2f Hz\n', pitch_E(frame_num));
fprintf('Pitch period at E: %.2f samples\n', ...
    pitch_period_E(frame_num));


%% Task 4 - Compare Pitch Estimates

figure;

plot(pitch, 'DisplayName', 'A - Original');

hold on;
plot(pitch_filtered, 'DisplayName', 'B - Filtered');
plot(pitch_E, 'DisplayName', 'E - Regenerated');

grid on;
legend;

xlabel('Frame Number');
ylabel('Fundamental Frequency (Hz)');
title('Pitch Estimation Before and After Regeneration');
ylim([80 420]);


%% Task 4 - Plot HPS at E

figure;

plot(fundamental_filters, hps_E(frame_num, :));
grid on;

xlabel('Filter Number');
ylabel('HPS');

title(sprintf( ...
    'Regenerated Speech HPS - Frame %d', frame_num));