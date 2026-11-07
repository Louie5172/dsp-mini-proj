%% Parameters

n = 512;

Qz = 1000000; % set to inf for this stuff
Qp = 8;

fs = 16e3;
dx = 3.5 / n;

% indices for analogue parameters
dist_idx = 1;
fp_idx = 2;
fn_idx = 3;
bwp_idx = 4;
bwz_idx = 5;

% indices for digital parameters
theta_p_idx = 1;
theta_z_idx = 2;
dbwp_idx = 3;
dbwz_idx = 4;

% indices for coefficients
Gp_idx = 1;
b1_idx = 2;
b2_idx = 3;
Gz_idx = 4;
a1_idx = 5;
a2_idx = 6;

% initialise matrices
params = zeros(n, 5);
digi = zeros(n, 4);
coeff = zeros(n, 6);


%% Frequency Calculations

filter_no = 1:n;

arr = 0:(n - 1);
params(:, dist_idx) = arr' * dx;

params(:, fp_idx) = 7500 * 10.^(-0.667 * params(:, dist_idx));

params(1, fn_idx) = 7580;
for i = 2:n
    params(i, fn_idx) = params(i - 1, fp_idx) * 0.985;
end


%% Bandwidth Calculations

for i = 1:n
    params(i, bwp_idx) = params(i, fp_idx) / Qp;
    params(i, bwz_idx) = params(i, fn_idx) / Qz;

    digi(i, theta_p_idx) = params(i, fp_idx) * 2 * pi / fs;
    digi(i, theta_z_idx) = params(i, fn_idx) * 2 * pi / fs;

    digi(i, dbwp_idx) = params(i, bwp_idx) * 2 * pi / fs;
    digi(i, dbwz_idx) = params(i, bwz_idx) * 2 * pi / fs;
end


%% Coefficient Calculations

for i = 1:n
    r = 1 - digi(i, dbwp_idx) / 2;

    coeff(i, b1_idx) = 2 * r * cos(digi(i, theta_p_idx));
    coeff(i, b2_idx) = r^2;

    rz = 1 - digi(i, dbwz_idx) / 2;
    coeff(i, a1_idx) = 2 * rz * cos(digi(i, theta_z_idx));
    coeff(i, a2_idx) = rz^2;

    coeff(i, Gp_idx) = 1 - coeff(i, b1_idx) + coeff(i, b2_idx);
    coeff(i, Gz_idx) = 1 / (1 - coeff(i, a1_idx) + coeff(i, a2_idx));
end


%% Impulse Simulation

N = fs;

x = [1; zeros(N-1, 1)];

displacement = zeros(N, n);
pressure = zeros(N, n);

y = x;

for i = 1:n
    b = coeff(i, Gz_idx) * [1, -coeff(i, a1_idx), coeff(i, a2_idx)];
    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement(:, i) = y;

    y = filter(b, 1, y);
    pressure(:, i) = y;
end


%% Spatial Differentiation - Impulse

spatial_diff1 = displacement(:, 1:end-1) - displacement(:, 2:end);
spatial_diff2 = spatial_diff1(:, 1:end-1) - spatial_diff1(:, 2:end);


%% FFT - Impulse Response

Y_disp = fft(displacement);
Y_disp1 = fft(spatial_diff1);
Y_disp2 = fft(spatial_diff2);

Y_press = fft(pressure);


%% Plot Impulse Responses

sections = [88, 226, 422];
samples = [300, 1500, 5000];

matrices = {displacement, spatial_diff1, spatial_diff2};
names = {'Displacement', 'First Difference', 'Second Difference'};

figure;

for i = 1:3
    for j = 1:3
        subplot(3, 3, (i-1)*3 + j);

        matrix = matrices{j};
        plot(1:samples(i), matrix(1:samples(i), sections(i)));

        title(sprintf('Filter %d - %s', sections(i), names{j}));
        xlabel('Sample Number');
        ylabel('Amplitude');
    end
end


%% Plot FFT Responses

section = 192;

matrices = {Y_disp, Y_disp1, Y_disp2};
names = {'No Spatial Difference', 'First Difference', 'Second Difference'};

f = (0:floor(N/2)) * fs / N;

figure;

for i = 1:3
    subplot(3, 1, i);

    matrix = matrices{i};

    mag = abs(matrix(1:length(f), section));

    % Normalise so the maximum is 0 dB
    mag = mag / max(mag);

    mag_db = 20 * log10(max(mag, eps));

    semilogx(f(2:end), mag_db(2:end));

    xlim([1, 1000]);
    ylim([-100, 0]);
    grid on;

    title(sprintf('Filter %d - %s', section, names{i}));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
end


%% Sinusoidal Simulation

f_sin = 1000;

t = (0:N-1)' / fs;
x_sin = sin(2*pi*f_sin*t);

displacement_sin = zeros(N, n);
pressure_sin = zeros(N, n);

y = x_sin;

for i = 1:n
    b = coeff(i, Gz_idx) * [1, -coeff(i, a1_idx), coeff(i, a2_idx)];
    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement_sin(:, i) = y;

    y = filter(b, 1, y);
    pressure_sin(:, i) = y;
end

spatial_sin1 = displacement_sin(:, 1:end-1) - displacement_sin(:, 2:end);
spatial_sin2 = spatial_sin1(:, 1:end-1) - spatial_sin1(:, 2:end);


%% Plot Sinusoidal Spatial Response

sample = 1000;

figure;

subplot(2, 1, 1);
plot(1:size(displacement_sin, 2), displacement_sin(sample, :));
xlim([0, 500]);
xlabel('Filter Number');
ylabel('Displacement');
title('Without Spatial Differentiation');

subplot(2, 1, 2);
plot(1:size(spatial_sin2, 2), spatial_sin2(sample, :));
xlim([0, 500]);
xlabel('Filter Number');
ylabel('Displacement');
title('After Two Spatial Differentiations');


%% Plot Pressure and Spatial Difference Envelopes

matrices = {pressure_sin, spatial_sin2};
num_samples_list = [8, 4];

for j = 1:2
    matrix = matrices{j};
    num_samples = num_samples_list(j);

    start_sample = round(0.5 * fs);
    samples = start_sample + (16 / num_samples) * (0:(num_samples - 1));

    figure;
    hold on;

    for i = samples
        plot(1:size(matrix, 2), matrix(i, :));
    end

    hold off;

    xlabel('Filter Number');
    ylabel('Amplitude');
    title('Envelope - 1 kHz Input');
end


%% Final Validation - Three-Tone Signal

N_test = 2000;
t_test = (0:N_test-1)' / fs;

% Three equal-amplitude sinusoids
f1 = 200;
f2 = 1000;
f3 = 4000;

x_test = sin(2*pi*f1*t_test) ...
       + sin(2*pi*f2*t_test) ...
       + sin(2*pi*f3*t_test);

% Run through cochlear filter bank
displacement_test = zeros(N_test, n);
pressure_test = zeros(N_test, n);

y = x_test;

for i = 1:n

    b = coeff(i, Gz_idx) * ...
        [1, -coeff(i, a1_idx), coeff(i, a2_idx)];

    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement_test(:, i) = y;

    y = filter(b, 1, y);
    pressure_test(:, i) = y;
end

% Spatial differentiation
spatial_test1 = displacement_test(:, 1:end-1) ...
              - displacement_test(:, 2:end);

spatial_test2 = spatial_test1(:, 1:end-1) ...
              - spatial_test1(:, 2:end);

%% Final Validation - Spatial Displacement Plots

sample_test = 1800;

figure;

subplot(3,1,1);
plot(1:n, displacement_test(sample_test, :));
grid on;
xlabel('Section Number');
ylabel('Displacement');
title('Before Spatial Differentiation');

subplot(3,1,2);
plot(1:n-1, spatial_test1(sample_test, :));
grid on;
xlabel('Section Number');
ylabel('First Difference');
title('After One Stage of Spatial Differentiation');

subplot(3,1,3);
plot(1:n-2, spatial_test2(sample_test, :));
grid on;
xlabel('Section Number');
ylabel('Second Difference');
title('After Two Stages of Spatial Differentiation');

sgtitle(sprintf('Three-Tone Cochlear Response at Sample %d', ...
    sample_test));


%% Final Validation - Frequency Selectivity

% Select 10 sections across the cochlear model
sections = round(linspace(20, 490, 10));

% FFT of impulse responses along the time dimension
N_fft = size(displacement, 1);

Y_disp = fft(displacement, N_fft, 1);
Y_disp2 = fft(spatial_diff2, N_fft, 1);

% Positive frequency axis
f = (0:floor(N_fft/2)) * fs/N_fft;

figure;

%% Before Spatial Differentiation

subplot(2,1,1);
hold on;

for i = 1:length(sections)

    section = sections(i);

    mag = abs(Y_disp(1:length(f), section));

    % Normalise each response to 0 dB
    mag = mag / max(mag);

    mag_db = 20*log10(max(mag, eps));

    semilogx(f(2:end), mag_db(2:end), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('Section %d', section));
end

grid on;
xlim([10 7500]);
ylim([-100 5]);

xlabel('Frequency (Hz)');
ylabel('Normalised magnitude (dB)');
title('Before Spatial Differentiation');
legend('Location', 'bestoutside');


%% After Two Stages of Spatial Differentiation

subplot(2,1,2);
hold on;

for i = 1:length(sections)

    section = sections(i);

    mag = abs(Y_disp2(1:length(f), section));

    % Normalise each response to 0 dB
    mag = mag / max(mag);

    mag_db = 20*log10(max(mag, eps));

    semilogx(f(2:end), mag_db(2:end), ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('Section %d', section));
end

grid on;
xlim([10 7500]);
ylim([-100 5]);

xlabel('Frequency (Hz)');
ylabel('Normalised magnitude (dB)');
title('After Two Stages of Spatial Differentiation');

sgtitle('Frequency Selectivity of Cochlear Model');
