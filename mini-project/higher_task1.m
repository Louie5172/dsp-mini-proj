%% Higher Level Task 1 - Generate Four Sinusoids

freqs = [60, 200, 1000, 2500];

num_samples = 1600;

t_input = (0:num_samples-1)' / fs;
x_input = zeros(num_samples, 1);

for f_sin = freqs
    x_input = x_input + sin(2*pi*f_sin*t_input);
end

% repeat the pattern for 1 second to allow settling
x_mix = repmat(x_input, ceil(fs/num_samples), 1);
x_mix = x_mix(1:fs);


%% Higher Level Task 1 - Filter Bank

N_mix = length(x_mix);

displacement_mix = zeros(N_mix, n);
pressure_mix = zeros(N_mix, n);

y = x_mix;

for i = 1:n
    b = coeff(i, Gz_idx) * [1, -coeff(i, a1_idx), coeff(i, a2_idx)];
    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement_mix(:, i) = y;

    y = filter(b, 1, y);
    pressure_mix(:, i) = y;
end

spatial_mix1 = displacement_mix(:, 1:end-1) - displacement_mix(:, 2:end);
spatial_mix2 = spatial_mix1(:, 1:end-1) - spatial_mix1(:, 2:end);


%% Higher Level Task 1 - Inner Hair Cell

fc = 30;
c_0 = exp(-2*pi*fc/fs);

% half-wave rectification
e_matrix = max(spatial_mix2, 0);

% first-order low-pass filter, zero initial conditions
energy = filter(1-c_0, [1, -c_0], e_matrix);

%% Higher Level Task 1 - Plot Analyser Output

% read after initial transients
sample = round(0.8 * fs) + 1;

energy_at_time = energy(sample, :);

energy_db = 20 * log10(max(energy_at_time, realmin));
energy_db = energy_db - max(energy_db);

figure;
plot(1:length(energy_db), energy_db);
grid on;

xlabel('Filter Number');
ylabel('Relative Output (dB)');
title('Inner Hair Cell Spectrum Analyser');

% identify strongest local maxima
peaks = find(energy_db(2:end-1) > energy_db(1:end-2) & ...
    energy_db(2:end-1) >= energy_db(3:end)) + 1;

[~, order] = sort(energy_db(peaks), 'descend');
peaks = peaks(order(1:min(4, length(order))));

hold on;
plot(peaks, energy_db(peaks), 'ro');
hold off;


%% Higher Level Task 1 - Input FFT

N_fft = length(x_input);

Y_input = fft(x_input);

f = (0:floor(N_fft/2)) * fs / N_fft;
mag = abs(Y_input(1:length(f))) / N_fft;
mag(2:end-1) = 2 * mag(2:end-1);

figure;
plot(f, mag);
grid on;

xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Input Signal FFT Magnitude');
xlim([0, 3000]);


%% Higher Level Task 1 - Output Sampling Comparison

section = 192;

energy_trace = energy(:, section);

start_sample = round(0.5 * fs) + 1;

samples16 = start_sample:round(0.016*fs):length(energy_trace);
samples8 = start_sample:round(0.008*fs):length(energy_trace);

figure;

subplot(2, 1, 1);
plot((samples16-1)/fs, energy_trace(samples16), '-o');
xlabel('Time (s)');
ylabel('Hair Cell Output');
title('16 ms Readout - 62.5 Hz');

subplot(2, 1, 2);
plot((samples8-1)/fs, energy_trace(samples8), '-o');
xlabel('Time (s)');
ylabel('Hair Cell Output');
title('8 ms Readout - 125 Hz');


%% Higher Level Task 1 - Frequency Resolution

freqs_wide = [1000, 1250];
freqs_close = [1000, 1050];

freq_sets = {freqs_wide, freqs_close};

db_results = cell(1, 2);

for test = 1:2

    freqs_test = freq_sets{test};

    % generate short input pattern
    x_test = zeros(num_samples, 1);

    for f_sin = freqs_test
        x_test = x_test + sin(2*pi*f_sin*t_input);
    end

    % repeat input
    x_test = repmat(x_test, ceil(fs/num_samples), 1);
    x_test = x_test(1:fs);

    % filter bank
    y = x_test;
    displacement_test = zeros(length(x_test), n);

    for i = 1:n
        b = coeff(i, Gz_idx) * [1, -coeff(i, a1_idx), coeff(i, a2_idx)];
        a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

        y = filter(coeff(i, Gp_idx), a, y);
        displacement_test(:, i) = y;

        y = filter(b, 1, y);
    end

    % second spatial difference
    spatial_test1 = displacement_test(:, 1:end-1) - displacement_test(:, 2:end);
    spatial_test2 = spatial_test1(:, 1:end-1) - spatial_test1(:, 2:end);

    % inner hair cell
    e_test = max(spatial_test2, 0);
    energy_test = filter(1-c_0, [1, -c_0], e_test);

    % convert to relative dB
    energy_at_time = energy_test(sample, :);
    db = 20 * log10(max(energy_at_time, realmin));

    db_results{test} = db - max(db);
end


%% Higher Level Task 1 - Plot Frequency Resolution

figure;

for test = 1:2
    subplot(2, 1, test);

    plot(1:length(db_results{test}), db_results{test});
    grid on;

    xlim([150, 220]);
    ylim([-60, 5]);

    xlabel('Filter Number');
    ylabel('Relative Output (dB)');

    freqs_test = freq_sets{test};
    title(sprintf('%d Hz and %d Hz', freqs_test(1), freqs_test(2)));
end
