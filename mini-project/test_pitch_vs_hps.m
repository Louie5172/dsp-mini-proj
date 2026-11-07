function test_pitch_vs_hps()
% TEST_PITCH_VS_HPS
% Compare MATLAB Audio Toolbox pitch estimates with the ELEC3104
% 512-section cochlear filter bank + inner hair cells + harmonic product.
%
% Put this file in the same folder as:
%   Female16K.wav
%   Harmonics_of_filters.csv
% Then run: test_pitch_vs_hps
%
% The function has its own workspace, so a variable called 'pitch' in your
% other MATLAB script will not shadow MATLAB's pitch() function.

%% Settings
useSynthetic = false;          % true: test a known 205 Hz, 4-harmonic input
knownF0 = 205;
frameNum = 20;                 % inspect the same frame as your assignment
Qp = 8;
Qz = 10000000;                      % try 16 or inf, as required by your assignment
fc = 30;
fsExpected = 16000;
n = 512;
dx = 0.00684;

%% Load speech (or construct a controlled test)
thisFolder = fileparts(mfilename('fullpath'));
if useSynthetic
    fs = fsExpected;
    t = (0:2*fs-1)'/fs;
    s = sin(2*pi*knownF0*t) ...
      + 0.8*sin(2*pi*2*knownF0*t) ...
      + 0.6*sin(2*pi*3*knownF0*t) ...
      + 0.4*sin(2*pi*4*knownF0*t);
    signalLabel = sprintf('Synthetic %.0f Hz harmonic signal', knownF0);
else
    audioPath = fullfile(thisFolder, 'Female16K.wav');
    assert(isfile(audioPath), 'Could not find Female16K.wav beside this script.');
    [s, fs] = audioread(audioPath);
    if size(s,2) > 1
        s = mean(s,2);  % mono
    end
    signalLabel = 'Female16K.wav';
end
assert(fs == fsExpected, 'Expected a 16000 Hz WAV file; actual sample rate is %g Hz.', fs);
s = s(:);
frameSize = 320;               % 20 ms at 16 kHz
numFrames = floor(numel(s)/frameSize);
assert(frameNum >= 1 && frameNum <= numFrames, 'Requested frame is outside the recording.');
frameEnds = (1:numFrames)'*frameSize;
frameTime = frameEnds/fs;       % right edge of each 20 ms frame

%% Built-in pitch estimators (Audio Toolbox)
assert(exist('pitch', 'file') == 2, ...
    'MATLAB pitch() was not found. Check that Audio Toolbox is installed.');

% Exact 20 ms, non-overlapping frames: aligns directly with your frame 20.
[f0NCF20, loc20] = pitch(s, fs, 'Method', 'NCF', ...
    'Range', [90 400], 'WindowLength', 320, 'OverlapLength', 0);

% Longer analysis windows are useful as a less noisy reference.
% Their returned locations are the LAST sample in each window.
[f0NCF40, loc40] = pitch(s, fs, 'Method', 'NCF', ...
    'Range', [90 400], 'WindowLength', 640, 'OverlapLength', 320);
[f0LHS40, locLHS] = pitch(s, fs, 'Method', 'LHS', ...
    'Range', [90 400], 'WindowLength', 640, 'OverlapLength', 320);

[~, j40] = min(abs(loc40 - frameEnds(frameNum)));
[~, jLHS] = min(abs(locLHS - frameEnds(frameNum)));

fprintf('\n--- %s; frame %d (%.3f-%.3f s) ---\n', ...
    signalLabel, frameNum, (frameNum-1)*frameSize/fs, frameNum*frameSize/fs);
fprintf('Built-in NCF, 20 ms: %.2f Hz\n', f0NCF20(frameNum));
fprintf('Built-in NCF, 40 ms: %.2f Hz\n', f0NCF40(j40));
fprintf('Built-in LHS, 40 ms: %.2f Hz\n', f0LHS40(jLHS));

%% Rebuild your exact 512-section cascade coefficients
fp = 7500 * 10.^(-0.667*(0:n-1)'*dx);
fn = [7580; 0.985*fp(1:end-1)];
bwp = fp/Qp;
bwz = fn/Qz;                 % = 0 if Qz = inf

rp = 1 - pi*bwp/fs;
rz = 1 - pi*bwz/fs;
b1 = 2*rp.*cos(2*pi*fp/fs);
b2 = rp.^2;
a1 = 2*rz.*cos(2*pi*fn/fs);
a2 = rz.^2;
Gp = 1 - b1 + b2;
Gz = 1./(1 - a1 + a2);

%% Cochlear filter bank: record displacement BEFORE each zero section
N = numel(s);
displacement = zeros(N,n);
y = s;
for i = 1:n
    y = filter(Gp(i), [1, -b1(i), b2(i)], y);
    displacement(:,i) = y;
    y = filter(Gz(i)*[1, -a1(i), a2(i)], 1, y);
end

%% Second spatial difference, rectification and 30 Hz hair-cell LPF
spatial2 = displacement(:,1:end-2) - 2*displacement(:,2:end-1) ...
         + displacement(:,3:end);
rectified = max(spatial2,0);
c0 = exp(-2*pi*fc/fs);
energy = filter(1-c0, [1, -c0], rectified);
energyFrames = energy(frameEnds,:);

%% HPS using the SAME harmonic filter lookup as your assignment
csvPath = fullfile(thisFolder,'Harmonics_of_filters.csv');
assert(isfile(csvPath), 'Could not find Harmonics_of_filters.csv beside this script.');
harmonics = readmatrix(csvPath);
assert(size(harmonics,2) >= 5, 'Harmonics CSV must contain at least five columns.');
harmonics = harmonics(all(isfinite(harmonics(:,1:5)),2),:);
fundamentalFreqs = harmonics(:,1);
fundamentalFilters = harmonics(:,2);
harmonicFilters = harmonics(:,2:5);
assert(all(harmonicFilters(:) == round(harmonicFilters(:))) && ...
       all(harmonicFilters(:) >= 1) && all(harmonicFilters(:) <= size(energy,2)), ...
       'Harmonic CSV references filter indices outside the 510 hair-cell channels.');

hpsAll = zeros(numFrames, size(harmonicFilters,1));
f0HPS = nan(numFrames,1);
for j = 1:numFrames
    v = energyFrames(j,:);
    hpsAll(j,:) = prod(v(harmonicFilters),2)';
    [~, idx] = max(hpsAll(j,:));
    f0HPS(j) = fundamentalFreqs(idx);
end
fprintf('Cochlear HPS, 20 ms: %.2f Hz\n',f0HPS(frameNum));

% Show the four actual channel energies for candidates near 205 and 220 Hz.
v = energyFrames(frameNum,:);
for candidate = [205,220]
    [~, idx] = min(abs(fundamentalFreqs-candidate));
    filters = harmonicFilters(idx,:);
    values = v(filters);
    fprintf('Candidate %.1f Hz: filters [%s], energies [%s], HPS = %.6g\n', ...
        fundamentalFreqs(idx), num2str(filters), num2str(values,'%.4g '), prod(values));
end

%% Figure 1: waveform and pitch estimates
figure('Name','Built-in MATLAB pitch vs cochlear HPS');
subplot(2,1,1);
plot((0:N-1)/fs,s);
xline(frameNum*frameSize/fs,'--','Frame 20 end');
grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title(signalLabel);

subplot(2,1,2);
plot(loc20/fs,f0NCF20,'-o','DisplayName','MATLAB NCF (20 ms)'); hold on;
plot(loc40/fs,f0NCF40,'-','LineWidth',1.5,'DisplayName','MATLAB NCF (40 ms)');
plot(locLHS/fs,f0LHS40,'-','LineWidth',1.5,'DisplayName','MATLAB LHS (40 ms)');
plot(frameTime,f0HPS,'-o','DisplayName','Cochlear HPS (20 ms)');
yline(knownF0,':','205 Hz guide');
xline(frameNum*frameSize/fs,'--');
hold off; grid on;
xlabel('Time (s)'); ylabel('Pitch (Hz)');
title('Pitch estimates (compare around frame 20)'); ylim([80 420]); legend('Location','best');

%% Figure 2: diagnostic plot of spatial response, energy and HPS at frame 20
figure('Name','Frame 20 cochlear stage diagnostics');
subplot(3,1,1);
plot(1:size(spatial2,2),spatial2(frameEnds(frameNum),:));
xlim([280 390]); grid on;
xlabel('Filter number'); ylabel('Displacement');
title(sprintf('Second spatial difference, frame %d end',frameNum));

subplot(3,1,2);
plot(1:size(energy,2),energyFrames(frameNum,:));
xlim([280 390]); grid on;
xlabel('Filter number'); ylabel('Hair-cell output');
title('Hair-cell output');

subplot(3,1,3);
plot(fundamentalFreqs,hpsAll(frameNum,:));
xline(205,'--','205 Hz'); xline(220,':','220 Hz');
grid on; xlabel('Candidate fundamental (Hz)'); ylabel('HPS');
title(sprintf('HPS: peak at %.1f Hz',f0HPS(frameNum)));
end
