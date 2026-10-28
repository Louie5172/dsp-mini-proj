% will be used for time domain simulation:

% include other script
find_param

N = fs; % 1 second of samples

% impulse input
x = [1; zeros(N-1, 1)];

% simulate impulse
[displacement, pressure] = simulate_filters(x, coeff, n, ...
    Gp_idx, b1_idx, b2_idx, Gz_idx, a1_idx, a2_idx);

% spatial differentiation
spatial_diff1 = displacement(:, 1:end-1) - displacement(:, 2:end);
spatial_diff2 = spatial_diff1(:, 1:end-1) - spatial_diff1(:, 2:end);

% FFT of impulse responses
Y_disp = fft(displacement);
Y_disp1 = fft(spatial_diff1);
Y_disp2 = fft(spatial_diff2);

Y_press = fft(pressure);


% sinusoidal input
f_sin = 1000;
t = (0:N-1)' / fs;
x_sin = sin(2*pi*f_sin*t);

% simulate sinusoid
[displacement_sin, pressure_sin] = simulate_filters(x_sin, coeff, n, Gp_idx, b1_idx, b2_idx, Gz_idx, a1_idx, a2_idx);

% spatial differentiation for sinusoid
spatial_sin1 = displacement_sin(:, 1:end-1) - displacement_sin(:, 2:end);
spatial_sin2 = spatial_sin1(:, 1:end-1) - spatial_sin1(:, 2:end);


% plotting

% impulse response
% plot_res(displacement, 88, 500);
% plot_res(spatial_diff2, 226, 5000);

% FFT magnitude responses
% plot_fft(Y_disp, 192, fs);
% plot_fft(Y_disp1, 192, fs);
% plot_fft(Y_disp2, 192, fs);

% sinusoidal response at section 192
% plot_res(displacement_sin, 192, 1000);

% sinusoidal displacement along BM at sample 1000
plot_res(displacement_sin, 222, 1000);
% plot_spatial(spatial_sin2, 1000, dx);


%% functions

% simulate cascaded filters
function [displacement, pressure] = simulate_filters(x, coeff, n, Gp_idx, b1_idx, b2_idx, Gz_idx, a1_idx, a2_idx)

    N = length(x);
    y = x;

    displacement = zeros(N, n);
    pressure = zeros(N, n);

    for i = 1:n
        b = coeff(i, Gz_idx) * [1, -coeff(i, a1_idx), coeff(i, a2_idx)];
        a = [1, -coeff(i, b1_idx), coeff(i, b2_idx)];

        % displacement
        y = filter(coeff(i, Gp_idx), a, y);
        displacement(:, i) = y;

        % pressure
        y = filter(b, 1, y);
        pressure(:, i) = y;
    end
end


% plot impulse response
function plot_res(matrix, section, samples)
    figure;
    plot(1:samples, matrix(1:samples, section));
    xlabel('Sample Number');
    ylabel('Amplitude');
end


% plot FFT magnitude response
function plot_fft(matrix, section, fs)
    N = size(matrix, 1);

    f = (0:floor(N/2)) * fs / N;
    mag = abs(matrix(1:length(f), section));

    figure;
    semilogx(f(2:end), 20*log10(max(mag(2:end), eps)));

    xlim([10 fs/2]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
end


% convert frequency to filter number
function num = freq_to_num(freq, dx)
    x = log10(freq / 7500) / -0.667;
    num = round(x / dx) + 1;
end