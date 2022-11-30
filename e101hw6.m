%% E101 HW6 MATLAB code
% Rafael Burger (rburger@hmc.edu) and Raja Batra (rbatra@hmc.edu)
% 11/29/2022

load("RVCradle.mat"); 

fs = 5000;    % sampling frequency of data in RVCradle
order = 3000; % order for non-Fourier PSD estimation
N = 800000;    % section of data to analyze for non-Fourier PSD estimation

blocktime = 10;               % block time in seconds (to average FRF)
windowlen = blocktime * fs;   % convert block time from seconds to sample num
overlap = 0.5 * windowlen;    % block overlap for use in calculating FRFs
w = hann(windowlen);          % hanning window of specified windowlen for use calculating FRFs

L2= L2Flow2S(1:N,2);       % truncate data to allow faster PSD estimation
L1 = L1Flow2S(1:N,2);
C = CFlow2S(1:N,2);
R1 = R1Flow2S(1:N,2);
R2 = R2Flow2S(1:N,2);

[psdL2, fL2] = pburg(L2, order, N, fs);   % perform PSD estimation
"completed L2"
[psdL1, fL1] = pburg(L1, order, N, fs);
"completed L1"
[psdC, fC]   = pburg(C, order, N, fs);
"completed C"
[psdR1, fR1] = pburg(R1, order, N, fs);
"completed R1"
[psdR2, fR2] = pburg(R2, order, N, fs);
"completed R2"

figure(1);
hold on;
plot(fL2, log10(psdL2));
plot(fL1, log10(psdL1));
plot(fC,  log10(psdC));
plot(fR1, log10(psdR1));
plot(fR2, log10(psdR2));
legend("L2", "L1", "C", "R1", "R2")
xlim([0 15]);


% [H_l1l2, f_l1l2] = tfestimate(L1Flow2S(:,2), L2Flow2S(:,2), w, overlap, windowlen, fs); % calculate FRFs
% [H_l1c, f_l1c]   = tfestimate(L1Flow2S(:,2), CFlow2S(:,2), w, overlap, windowlen, fs);
% [H_l1r1, f_l1r1] = tfestimate(L1Flow2S(:,2), R1Flow2S(:,2), w, overlap, windowlen, fs);
% [H_l1r2, f_l1r2] = tfestimate(L1Flow2S(:,2), R2Flow2S(:,2), w, overlap, windowlen, fs);
% [H_l1l1, f_l1l1] = tfestimate(L1Flow2S(:,2), L1Flow2S(:,2), w, overlap, windowlen, fs);

    

