
%% Debug 205 Hz - Full Cochlear Cascade

clear; clc; close all;

%% Parameters

n = 512;
Qp = 8;
Qz = 16;
fs = 16000;
dx = 0.00684;
f_test = 205;

% Parameter matrix indices
dist_idx = 1;
fp_idx = 2;
fn_idx = 3;
bwp_idx = 4;
bwz_idx = 5;

% Digital parameter indices
theta_p_idx = 1;
theta_z_idx = 2;
dbwp_idx = 3;
dbwz_idx = 4;

% Coefficient indices
Gp_idx = 1;
b1_idx = 2;
b2_idx = 3;
Gz_idx = 4;
a1_idx = 5;
a2_idx = 6;

params = zeros(n, 5);
digi = zeros(n, 4);
coeff = zeros(n, 6);


%% Frequency Calculations

params(:, dist_idx) = (0:n-1)' * dx;
params(:, fp_idx) = 7500 * 10.^(-0.667 * params(:, dist_idx));

params(1, fn_idx) = 7580;

for i = 2:n
    params(i, fn_idx) = params(i-1, fp_idx) * 0.985;
end


%% Bandwidth Calculations

for i = 1:n
    params(i, bwp_idx) = params(i, fp_idx) / Qp;
    params(i, bwz_idx) = params(i, fn_idx) / Qz;

    digi(i, theta_p_idx) = 2*pi*params(i, fp_idx)/fs;
    digi(i, theta_z_idx) = 2*pi*params(i, fn_idx)/fs;

    digi(i, dbwp_idx) = 2*pi*params(i, bwp_idx)/fs;
    digi(i, dbwz_idx) = 2*pi*params(i, bwz_idx)/fs;
end


%% Coefficient Calculations

for i = 1:n
    rp = 1 - digi(i, dbwp_idx)/2;
    rz = 1 - digi(i, dbwz_idx)/2;

    coeff(i, b1_idx) = 2*rp*cos(digi(i, theta_p_idx));
    coeff(i, b2_idx) = rp^2;

    coeff(i, a1_idx) = 2*rz*cos(digi(i, theta_z_idx));
    coeff(i, a2_idx) = rz^2;

    coeff(i, Gp_idx) = 1 - coeff(i, b1_idx) + coeff(i, b2_idx);
    coeff(i, Gz_idx) = 1 / (1 - coeff(i, a1_idx) + coeff(i, a2_idx));
end


%% Steady-State Transfer Functions at 205 Hz

% Evaluate each section at z^-1 = exp(-j*omega)
z_inv = exp(-1j * 2*pi*f_test/fs);

Hp = zeros(1, n);
Hz = zeros(1, n);
Hdisp = zeros(1, n);
Hpressure = zeros(1, n);

H_previous = 1;

for i = 1:n

    Hp(i) = coeff(i, Gp_idx) / ...
        (1 - coeff(i, b1_idx)*z_inv + coeff(i, b2_idx)*z_inv^2);

    Hz(i) = coeff(i, Gz_idx) * ...
        (1 - coeff(i, a1_idx)*z_inv + coeff(i, a2_idx)*z_inv^2);

    % Displacement measured after pole stage
    Hdisp(i) = H_previous * Hp(i);

    % Pressure passed to next section
    H_previous = Hdisp(i) * Hz(i);
    Hpressure(i) = H_previous;
end


%% Gain Between Consecutive Sections

% |displacement(i+1)| / |displacement(i)|
gain_ratio = abs(Hz(1:end-1) .* Hp(2:end));

% Second spatial difference of complex frequency responses
Hdiff2 = Hdisp(1:end-2) - 2*Hdisp(2:end-1) + Hdisp(3:end);

% Find expected and actual peaks
[~, expected_filter] = min(abs(params(:, fp_idx) - f_test));
[~, disp_peak] = max(abs(Hdisp));
[~, diff_peak] = max(abs(Hdiff2));

fprintf('Input frequency: %.2f Hz\n', f_test);
fprintf('Qp = %g, Qz = %g\n', Qp, Qz);

fprintf('Nominal nearest channel: %d (%.2f Hz)\n', ...
    expected_filter, params(expected_filter, fp_idx));

fprintf('Displacement peak: %d (%.2f Hz nominal)\n', ...
    disp_peak, params(disp_peak, fp_idx));

fprintf('Second-difference peak: %d (%.2f Hz nominal)\n', ...
    diff_peak, params(diff_peak, fp_idx));


%% Figure 1 - Trace Where Response Begins to Decay

figure;

% displacement magnitude
subplot(3, 1, 1);
plot(1:n, abs(Hdisp), 'LineWidth', 1.2);
hold on;

xline(expected_filter, '--r', '205 Hz channel');
xline(disp_peak, '--k', 'Displacement peak');

hold off;
xlim([300 365]);
grid on;

ylabel('|D_i(205 Hz)|');
title('Steady-State Displacement Magnitude');


% gain between consecutive sections
subplot(3, 1, 2);
plot(1:n-1, gain_ratio, 'LineWidth', 1.2);
hold on;

yline(1, '--k', 'Unity gain');
xline(expected_filter, '--r');

hold off;
xlim([300 365]);
grid on;

ylabel('|D_{i+1}| / |D_i|');
title('Gain Between Consecutive Sections');


% second spatial difference
subplot(3, 1, 3);
plot(1:n-2, abs(Hdiff2), 'LineWidth', 1.2);
hold on;

xline(expected_filter, '--r', '205 Hz channel');
xline(diff_peak, '--k', 'Spatial peak');

hold off;
xlim([300 365]);
grid on;

xlabel('Filter Number');
ylabel('|Second Difference|');
title('Steady-State Second Spatial Difference');


%% Time-Domain Verification

% Generate pure 205 Hz sinusoid
N = fs;
sample_idx = (0:N-1)';

x = sin(2*pi*f_test*sample_idx/fs);

displacement = zeros(N, n);

y = x;

for i = 1:n

    b = coeff(i, Gz_idx) * ...
        [1, -coeff(i, a1_idx), coeff(i, a2_idx)];

    a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

    y = filter(coeff(i, Gp_idx), a, y);
    displacement(:, i) = y;

    y = filter(b, 1, y);
end


%% Two Spatial Differentiations

spatial_diff1 = displacement(:, 1:end-1) - displacement(:, 2:end);
spatial_diff2 = spatial_diff1(:, 1:end-1) - spatial_diff1(:, 2:end);


%% Inner Hair Cell Model

fc = 30;
c_0 = exp(-2*pi*fc/fs);

% Half-wave rectification
e_matrix = max(spatial_diff2, 0);

% Low-pass filtering
energy = filter(1-c_0, [1, -c_0], e_matrix);


%% Calculate RMS and Mean Responses

% Use settled samples to avoid initial transients
settled = round(0.5*fs):round(0.9*fs);

disp_rms = sqrt(mean(displacement(settled, :).^2, 1));

diff_rms = sqrt(mean(spatial_diff2(settled, :).^2, 1));

energy_mean = mean(energy(settled, :), 1);

[~, ihc_peak] = max(energy_mean);

fprintf('Mean hair-cell peak: %d (%.2f Hz nominal)\n', ...
    ihc_peak, params(ihc_peak, fp_idx));


%% Figure 2 - Time Domain Verification

figure;

% Displacement RMS
subplot(3, 1, 1);
plot(1:n, disp_rms);
hold on;

xline(expected_filter, '--r', '205 Hz channel');

hold off;
xlim([300 365]);
grid on;

ylabel('RMS');
title('Displacement RMS');


% Second spatial difference RMS
subplot(3, 1, 2);
plot(1:n-2, diff_rms);
hold on;

xline(expected_filter, '--r', '205 Hz channel');

hold off;
xlim([300 365]);
grid on;

ylabel('RMS');
title('Second Spatial Difference RMS');


% Inner hair cell output
subplot(3, 1, 3);
plot(1:n-2, energy_mean);
hold on;

xline(expected_filter, '--r', '205 Hz channel');
xline(ihc_peak, '--k', 'Hair-cell peak');

hold off;
xlim([300 365]);
grid on;

xlabel('Filter Number');
ylabel('Mean Energy');
title('Mean Inner Hair Cell Output');
