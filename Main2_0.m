clc; clear ; close all;

%% --------------------- Task 0 -----------------
fprintf('\n=== TASK 0: Audio File ===\n\n');
%Reading Sound File & Defining main parameters: 
[x, Fs] = audioread("sample_audio_file.wav");   % returns x: Audio signal samples , Fs: Sampling frequency (Found in file prop.)
N = length(x);
df= Fs/N;
ts=1/Fs;
t = (0:N-1) / Fs;                 % time axis [nts]
f = (-Fs/2 : df : (Fs/2)-df);     % -Fs/2 < f < Fs/2  

% --------------- Extras --------------
%{
figure;
plot(t, x);
xlabel("Time (sec)");
ylabel("Amplitude");
title("Time Domain Sample Audio File x[n]");   % the plot will have 2 colors as the audio file is stereo not mono -> 2 vectors
grid on; 
%}
% --------------- End of Extras -------


%As we saw the plot had 2 colors -> audio is stereo 
%Converting the audio to mono -> to handle it as 1 vector signal samples
x_mono = mean(x, 2);


X = fft(x_mono);              % compute FFT
X_shifted = fftshift(X)*ts;   % center 0 frequency

figure; 

subplot(2,1,1);
plot(t, x_mono);
xlabel("Time (sec)");
ylabel("Amplitude");
title("Time Domain Sample Audio File monotonic x[n]");
grid on;

subplot(2,1,2);
plot(f, abs(X_shifted));
xlabel("Frequency (Hz)");
ylabel("|X(f)|");
title("Magnitude Spectrum of the Sample Audio File monotonic X[f]");
grid on;

Energy_time = sum(abs(x_mono).^2)*ts

Energy_freq = sum(abs(X_shifted).^2)*df

%% --------------------- Task 1 -----------------
fprintf('\n=== TASK 1: Digital Echo System ===\n\n');

% Given parameters
b_k = [1, 0.9, 0.8, 0.7];
D = 1000;

% Build FIR numerator: H(z)=1 + 0.9 z^-1000 + 0.8 z^-2000 + 0.7 z^-3000
num_H = zeros(1, 3*D + 1);
num_H(1)         = b_k(1);
num_H(D+1)       = b_k(2);
num_H(2*D+1)     = b_k(3);
num_H(3*D+1)     = b_k(4);

den_H = 1;   % FIR filter denominator

fprintf("Generating plots for Echo System H(z)...\n");
Sec3Tasks_2_1(' - Echo System H(z)', num_H, den_H, num_H, den_H, f, Fs, N, 3500 ,true);

y1 = filter(num_H, den_H, x_mono);

% Compute MSE for y1
MSE_y1 = calc_mse(y1, x_mono);

% Normalize and save the echoed audio
y1_normalized = y1 / max(abs(y1));
audiowrite('audio_with_echo.wav', y1_normalized, Fs);


% G(z) = 1 / H(z)  →  num_G = 1, den_G = num_H
num_G = 1;
den_G = num_H;

fprintf("Generating plots for Equalizer System G(z)...\n");
Sec3Tasks_2_1(' - Equalizer System G(z)', num_G, den_G, num_G, den_G, f, Fs, N , 20000, true);

y2 = filter(num_G, den_G, y1);

% Compute MSE for y2
MSE_y2 = calc_mse(y2, x_mono);

% Normalize and save
y2_normalized = y2 / max(abs(y2));
audiowrite('audio_equalized.wav', y2_normalized, Fs);

fprintf('\n== RESULTS ==\n');
fprintf('MSE (y1 - echoed):     %.6e\n', MSE_y1);
fprintf('MSE (y2 - equalized):  %.6e\n', MSE_y2);

%% --------------- Extras - Comparison -------------
%{
plot_samples = min(10000, N);

figure;
subplot(3,1,1);
plot(t(1:plot_samples), x_mono(1:plot_samples));
xlabel('Time (sec)');
ylabel('Amplitude');
title('Original Signal x[n]');
grid on;

subplot(3,1,2);
plot(t(1:plot_samples), y1(1:plot_samples));
xlabel('Time (sec)');
ylabel('Amplitude');
title('Signal with Echo y1[n]');
grid on;

subplot(3,1,3);
plot(t(1:plot_samples), y2(1:plot_samples));
xlabel('Time (sec)');
ylabel('Amplitude');
title('Equalized Signal y2[n]');
grid on;
%}


%% --------------------- Task 2: Design of Digital IIR LPF -----------------
fprintf('\n=== TASK 2: Design of Digital IIR LPF ===\n\n');

% filter specs
fp = 3000;      % passband edge in Hz
fstop = 4000;   % stopband edge in Hz
Ap = 1;         % max passband ripple (dB)
As = 50;        % min stopband attenuation (dB)

% normalize frequencies
Wp = fp / (Fs/2);       % passband edge
Ws = fstop / (Fs/2);    % stopband edge

% Design all four filter types
fprintf('=== Designing Filters ===\n');

% BUTTERWORTH FILTER
[n_butter, Wn_butter] = buttord(Wp, Ws, Ap, As);          % returns min order & cuff off freq
[b_butter, a_butter] = butter(n_butter, Wn_butter,'low'); %returns filter coefficients

% CHEBYSHEV TYPE I FILTER
[n_cheby1, Wn_cheby1] = cheb1ord(Wp, Ws, Ap, As);
[b_cheby1, a_cheby1] = cheby1(n_cheby1, Ap, Wn_cheby1,'low');

% CHEBYSHEV TYPE II FILTER
[n_cheby2, Wn_cheby2] = cheb2ord(Wp, Ws, Ap, As);
[b_cheby2, a_cheby2] = cheby2(n_cheby2, As, Wn_cheby2,'low');

% ELLIPTIC FILTER
[n_ellip, Wn_ellip] = ellipord(Wp, Ws, Ap, As);
[b_ellip, a_ellip] = ellip(n_ellip, Ap, As, Wn_ellip,'low');


% Filter Order Comparison Table
fprintf('\n========================================\n');
fprintf('     FILTER ORDER COMPARISON\n');
fprintf('========================================\n');
fprintf('%-25s | Order\n', 'Filter Type');
fprintf('-------------------------|-------\n');
fprintf('%-25s | %d\n', 'Butterworth', n_butter);
fprintf('%-25s | %d\n', 'Chebyshev Type I', n_cheby1);
fprintf('%-25s | %d\n', 'Chebyshev Type II', n_cheby2);
fprintf('%-25s | %d\n', 'Elliptic', n_ellip);
fprintf('========================================\n\n');

% Apply filters to audio and compute metrics
fprintf('=== Applying Filters to Audio Signal ===\n\n');

% Store filter data
filters = {'butter', 'cheby1', 'cheby2', 'ellip'};
filter_names = {'Butterworth', 'Chebyshev Type I', 'Chebyshev Type II', 'Elliptic'};
filter_coefs = { % a 1*4 array with each cell having a struct b and a, b for num coeff. and a for den coeff.
    struct('b', b_butter, 'a', a_butter);
    struct('b', b_cheby1, 'a', a_cheby1);
    struct('b', b_cheby2, 'a', a_cheby2);
    struct('b', b_ellip, 'a', a_ellip)
};

% storage for results
MSE_results = zeros(4, 1);         %empty 4*1 array
Energy_loss_results = zeros(4, 1); %empty 4*1 array

for i = 1:4
    fprintf('Processing: %s...\n', filter_names{i});
    
    b = filter_coefs{i}.b;
    a = filter_coefs{i}.a;
    
    % high order filters cause instability due to their very large/small
    % coefficients
    % convert to second-order section for numerical stability
    [z, p, k] = tf2zpk(b, a);
    [sos, g] = zp2sos(z, p, k);
    
    % Generate plots using SOS (Section 3 requirements)
    Sec3Tasks_2_1([' - ' filter_names{i} ' LPF'], sos, 1, b, a, f, Fs, N, 600, false, g); %pass num as SOS and den is ignored for second order sections
    
    % apply filter to audio
    y_filtered = filter(b, a, x_mono);
    
    % Compute MSE
    MSE_results(i) = calc_mse(y_filtered, x_mono);
    
    % Compute energy loss percentage
    Energy_loss_results(i) = calc_energylost(x_mono, y_filtered);
    
end

% Results Summary
fprintf('\n========================================\n');
fprintf('     FILTER PERFORMANCE COMPARISON\n');
fprintf('========================================\n');
fprintf('%-25s | MSE           | Energy Loss (%%)\n', 'Filter Type');
fprintf('-------------------------|---------------|----------------\n');
for i = 1:4
    fprintf('%-25s | %.6e | %.4f\n', filter_names{i}, MSE_results(i), Energy_loss_results(i));
end
fprintf('========================================\n\n');

% find minimum mse and energy lost
[min_MSE, idx_MSE] = min(MSE_results);
[min_loss, idx_loss] = min(Energy_loss_results);

fprintf('=== BEST FILTER ANALYSIS ===\n');
fprintf('Best filter (minimum MSE): %s\n', filter_names{idx_MSE});
fprintf('   MSE = %.6e\n\n', min_MSE);
fprintf('Best filter (minimum energy loss): %s\n', filter_names{idx_loss});
fprintf('   Energy Loss = %.4f%%\n\n', min_loss);

%%
%% --------------------- Task 3: Frequency Transformation using Pole-Zero Rotation -----------------
fprintf('\n=== TASK 3: Frequency Transformation using Pole-Zero Rotation ===\n\n');

%% Part A: Transform Butterworth LPF to HPF (rotate by π)
fprintf('=== Part A: Highpass Filter (HPF) via π rotation ===\n');

[z_lpf, p_lpf, k_lpf] = tf2zp(b_butter, a_butter);

% Rotate by π: multiply by e^(jπ) = -1
z_hpf = -z_lpf;
p_hpf = -p_lpf;
k_hpf = k_lpf;

[b_hpf, a_hpf] = zp2tf(z_hpf, p_hpf, k_hpf);

fprintf('Generating plots for HPF...\n');
[z, p, k] = tf2zpk(b_hpf, a_hpf);
[sos, g] = zp2sos(z, p, k);
Sec3Tasks_2_1(' - Butterworth HPF (π rotation)', sos, 1, b_hpf, a_hpf, f, Fs, N, 600, false, g);

%% Part B: Transform Butterworth LPF to BPF (rotate by π/2)
fprintf('\n=== Part B: Bandpass Filter (BPF) centered at π/2 ===\n');

% Manual pole-zero rotation for digital BPF
z_bpf_pos = z_lpf.*j;
z_bpf_neg = z_lpf./j;

p_bpf_pos = p_lpf.*j;     % Rotate by +π/2
p_bpf_neg = p_lpf./j;     % Rotate by -π/2

% Create two separate filters
z_bpf_1 = z_bpf_pos ;
p_bpf_1 = p_bpf_pos ;
k_bpf_1 = k_lpf;

z_bpf_2 = z_bpf_neg;
p_bpf_2 = p_bpf_neg;
k_bpf_2 = k_lpf;

b_bpf_1 = k_bpf_1*poly(z_bpf_1);
a_bpf_1 = poly(p_bpf_1);
b_bpf_2 = k_bpf_2*poly(z_bpf_2);
a_bpf_2 = poly(p_bpf_2);

% Parallel combination: H(z) = H1(z) + H2(z)
% This means: (b1/a1) + (b2/a2) = (b1*a2 + b2*a1) / (a1*a2)
% Parallel combination: H(z) = H1(z) + H2(z)
% (b1/a1) + (b2/a2) = (b1*a2 + a1*b2) / (a1*a2)
b_bpf = conv(b_bpf_1, a_bpf_2) + conv(a_bpf_1, b_bpf_2);
a_bpf = conv(a_bpf_1, a_bpf_2);

b_bpf = real(b_bpf);
a_bpf = real(a_bpf);

[z, p, k] = tf2zpk(b_bpf, a_bpf); 
[sos, g] = zp2sos(z,p,k);
fprintf('Generating plots for BPF...\n');
Sec3Tasks_2_1(' - Butterworth BPF (π/2 rotation)', sos, 1, b_bpf, a_bpf, f, Fs, N, 600, false, g);

%% ----------EXTRAS Comparison Plot: LPF vs HPF vs BPF------------------
fprintf('\n=== Generating comparison plots ===\n');

figure('Name', 'Task 3: Filter Comparison (LPF, HPF, BPF)', 'Position', [150, 150, 1400, 600]);

[z, p, k] = tf2zpk(b_butter, a_butter); 
[sos_lpf, g_lpf] = zp2sos(z,p,k);

[z, p, k] = tf2zpk(b_hpf, a_hpf); 
[sos_hpf, g_hpf] = zp2sos(z,p,k);

[z, p, k] = tf2zpk(b_bpf, a_bpf); 
[sos_bpf, g_bpf] = zp2sos(z,p,k);

H_lpf_t3 = g_lpf * freqz(sos_lpf, N, Fs, 'whole');
H_hpf_t3 = g_hpf*freqz(sos_hpf, N, Fs, 'whole');
H_bpf_t3 = g_bpf*freqz(sos_bpf, N, Fs, 'whole');

H_lpf_shift = fftshift(H_lpf_t3);
H_hpf_shift = fftshift(H_hpf_t3);
H_bpf_shift = fftshift(H_bpf_t3);

subplot(1,2,1);
hold on;
plot(f, 20*log10(abs(H_lpf_shift) + eps), 'b-', 'LineWidth', 2);
plot(f, 20*log10(abs(H_hpf_shift) + eps), 'r-', 'LineWidth', 2);
plot(f, 20*log10(abs(H_bpf_shift) + eps), 'g-', 'LineWidth', 2);
hold off;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response Comparison (Full View)');
legend('LPF (Original)', 'HPF (π rotation)', 'BPF (π/2 rotation)', 'Location', 'best');
grid on;
xlim([0, Fs/2]);
ylim([-100, 5]);

subplot(1,2,2);
hold on;
plot(f, 20*log10(abs(H_lpf_shift) + eps), 'b-', 'LineWidth', 2);
plot(f, 20*log10(abs(H_hpf_shift) + eps), 'r-', 'LineWidth', 2);
plot(f, 20*log10(abs(H_bpf_shift) + eps), 'g-', 'LineWidth', 2);
xline(3000, '--k', 'LPF fp', 'LineWidth', 1);
xline(4000, '--k', 'LPF fs', 'LineWidth', 1);
xline(Fs/4, '--m', 'BPF center', 'LineWidth', 1);
hold off;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Magnitude Response Comparison (Zoomed)');
legend('LPF (Original)', 'HPF (π rotation)', 'BPF (π/2 rotation)', 'Location', 'best');
grid on;
xlim([0, 15000]);
ylim([-80, 5]);

fprintf('\n=== End Of Project ===\n');
