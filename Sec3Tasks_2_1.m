function Sec3Tasks_2_1(String_title, num, den, num_zplane, den_zplane, Freq_Vector, F_sampling, N_Samples, imp_len ,show_zoom)
    % show_zoom: 
    %   - true: shows zoomed plots
    %   - false: shows only full-range plots
    
    %%---------------------------------- 1: Pole-Zero Plot (Separate Window) ----------------------------------
    figure('Name', ['Pole-Zero Plot ' String_title], 'NumberTitle','off');
    zplane(num_zplane, den_zplane);
    title(['Pole-Zero Plot ' String_title]);
    grid on;
    
    %%---------------------------------- Create Main Tabbed Figure ----------------------------------
    f = figure('Name', ['Filter Analysis ' String_title], 'NumberTitle', 'off', ...
               'Position',[100 50 1400 900]);
    tg = uitabgroup(f);
    tab2 = uitab(tg,'Title','Magnitude');
    tab3 = uitab(tg,'Title','Phase');
    tab4 = uitab(tg,'Title','Group Delay');
    tab5 = uitab(tg,'Title','Impulse Response');
    
    %%---------------------------------- Compute Frequency Data ----------------------------------
    if size(num, 2) == 6 %detects if num is in SOS format, which is a matrix with 6 columns
        % SOS format - only needs 2 arguments
        H = freqz(num, N_Samples, 'whole');                  % without whole: 0->fs/2 , with whole 0->fs
    else
        % Normal (b,a) format
        H = freqz(num, den, N_Samples, F_sampling, 'whole'); % without whole: 0->fs/2 , with whole 0->fs
    end
    H_shift = fftshift(H);
    x_mag_dB = 20 * log10(abs(H_shift) + eps);          % EPS is very small constant to avoid log 0
    x_phase  = angle(H_shift);
    phase_unwrapped = unwrap(x_phase);                  % removes discontinuities in phase at +-Pi
    
    
    % Group Delay
    % Computing group delay using logarithmic identity
   
    w = linspace(-pi, pi, N_Samples).';   %define freq vector from -pi to pi

    Hc = H_shift(:);                      %transforms vector into column vector to match

    
    dH_dw = gradient(Hc, w);              %computes derivative of H(e^jw) with respect to w

    % group delay
    gd = -imag(dH_dw ./ Hc);              % H'/H is the derivative of ln(H(e^jw)), taking the imaginary part of it represents derivative of phase

    % masks zeros and deep stopband 
    mag_thresh = 1e-6;
    gd(abs(Hc) < mag_thresh) = NaN;
    
    %%---------------------------------- 2: Magnitude TAB ----------------------------------
    p2 = uipanel('Parent',tab2);
    
    if show_zoom
        % Two subplots: full range + zoomed
        ax1 = subplot(2,1,1,'Parent',p2);
        plot(ax1, Freq_Vector, x_mag_dB);
        xlabel(ax1,'Frequency (Hz)');
        ylabel(ax1,'Magnitude (dB)');
        title(ax1,['Magnitude Response ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-F_sampling/2, F_sampling/2]);
        
        ax2 = subplot(2,1,2,'Parent',p2);
        plot(ax2, Freq_Vector, x_mag_dB);
        xlabel(ax2,'Frequency (Hz)');
        ylabel(ax2,'Magnitude (dB)');
        title(ax2,['Magnitude Response (Zoomed) ' String_title]);
        grid(ax2,'on');
        xlim(ax2,[-100,100]);
    else
        % Single plot: full range only
        ax1 = axes('Parent',p2);
        plot(ax1, Freq_Vector, x_mag_dB);
        xlabel(ax1,'Frequency (Hz)');
        ylabel(ax1,'Magnitude (dB)');
        title(ax1,['Magnitude Response ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-F_sampling/2, F_sampling/2]);
    end
    
    %%---------------------------------- 3: Phase TAB ----------------------------------
    p3 = uipanel('Parent',tab3);
    
    if show_zoom
        % Two subplots: full range + zoomed
        ax1 = subplot(2,1,1,'Parent',p3);
        plot(ax1, Freq_Vector, phase_unwrapped);
        xlabel(ax1,'Frequency (Hz)');
        ylabel(ax1,'Phase (rad)');
        title(ax1,['Phase Response ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-F_sampling/2, F_sampling/2]);
        
        ax2 = subplot(2,1,2,'Parent',p3);
        plot(ax2, Freq_Vector, phase_unwrapped);
        xlabel(ax2,'Frequency (Hz)');
        ylabel(ax2,'Phase (rad)');
        title(ax2,['Phase Response (Zoomed) ' String_title]);
        grid(ax2,'on');
        xlim(ax2,[-100,100]);
    else
        % Single plot: full range only
        ax1 = axes('Parent',p3);
        plot(ax1, Freq_Vector, phase_unwrapped);
        xlabel(ax1,'Frequency (Hz)');
        ylabel(ax1,'Phase (rad)');
        title(ax1,['Phase Response ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-F_sampling/2, F_sampling/2]);
    end
    
    %%---------------------------------- 4: Group Delay TAB ----------------------------------
    p4 = uipanel('Parent',tab4);
    
    if show_zoom
        % Two subplots: full range + zoomed
        ax1 = subplot(2,1,1,'Parent',p4);
        plot(ax1, w, gd);
        xlabel(ax1,'\omega (rad/sample)');
        ylabel(ax1,'Group Delay (samples)');
        title(ax1,['Group Delay ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-pi,pi]);
        
        ax2 = subplot(2,1,2,'Parent',p4);
        plot(ax2, w, gd);
        xlabel(ax2,'\omega (rad/sample)');
        ylabel(ax2,'Group Delay (samples)');
        title(ax2,['Group Delay (Zoomed) ' String_title]);
        grid(ax2,'on');
        xlim(ax2,[-0.01,0.01]);
    else
        % Single plot: full range only
        ax1 = axes('Parent',p4);
        plot(ax1, w, gd);
        xlabel(ax1,'\omega (rad/sample)');
        ylabel(ax1,'Group Delay (samples)');
        title(ax1,['Group Delay ' String_title]);
        grid(ax1,'on');
        xlim(ax1,[-pi,pi]);
    end
    
    %%---------------------------------- 5: Impulse Response TAB ----------------------------------
    if size(num, 2) == 6  % SOS format
        [h_imp, n_imp] = impz(num, imp_len);
    else  % Normal (b, a) format
        [h_imp, n_imp] = impz(num, den, imp_len);
    end
    
    p5 = uipanel('Parent',tab5);
    ax = axes('Parent',p5);
    stem(ax, n_imp, h_imp, 'filled', 'MarkerSize', 2);
    xlabel(ax,'n (samples)');
    ylabel(ax,'h[n]');
    title(ax,['Impulse Response ' String_title]);
    grid(ax,'on');
end