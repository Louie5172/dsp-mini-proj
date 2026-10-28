
% this is my dsp mini project

n = 512;

Qz = 16;
Qp = 8;

fs = 16e3;
dx = 0.00684;

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

% frequency calculations

filter_no = 1:n;

arr = 0:(n - 1);
params(:, dist_idx) = arr' * dx;

for i = 1:n
    params(i, fp_idx) = get_freq(params(i, dist_idx));
end

params(1, fn_idx) = 7580;
for i = 2:n
    params(i, fn_idx) = params(i - 1, fp_idx) * 0.985;
end

% bandwidth calculations

for i = 1:n
    params(i, bwp_idx) = params(i, fp_idx) / Qp;
    params(i, bwz_idx) = params(i, fn_idx) / Qz;

    % convert to the digital freq domain
    digi(i, theta_p_idx) = params(i, fp_idx) * 2 * pi / fs;
    digi(i, theta_z_idx) = params(i, fn_idx) * 2 * pi / fs;

    digi(i, dbwp_idx) = params(i, bwp_idx) * 2 * pi / fs;
    digi(i, dbwz_idx) = params(i, bwz_idx) * 2 * pi / fs;
end

% coefficient calculations

for i = 1:n
    r = 1 - digi(i, dbwp_idx) / 2;

    coeff(i, b1_idx) = 2 * r * cos(digi(i, theta_p_idx));
    coeff(i, b2_idx) = r^2;

    coeff(i, a1_idx) = 2 * cos(digi(i, theta_z_idx));
    coeff(i, a2_idx) = 1; % since we set r_z on the unit cicle so that it is always greater than r_p

    coeff(i, Gp_idx) = 1 - coeff(i, b1_idx) + coeff(i, b2_idx);
    coeff(i, Gz_idx) = 1 / (1 - coeff(i, a1_idx) + coeff(i, a2_idx));
end

% display
% format bank
% params
% digi
% coeff

function freq = get_freq(x)
    freq = 7500 * 10^(-0.667 * x);
end