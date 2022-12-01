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

% [psdL2, fL2] = pburg(L2, order, N, fs);   % perform PSD estimation
% "completed L2"
% [psdL1, fL1] = pburg(L1, order, N, fs);
% "completed L1"
% [psdC, fC]   = pburg(C, order, N, fs);
% "completed C"
% [psdR1, fR1] = pburg(R1, order, N, fs);
% "completed R1"
% [psdR2, fR2] = pburg(R2, order, N, fs);
% "completed R2"

% figure(1);
% hold on;
% plot(fL2, log10(psdL2));
% plot(fL1, log10(psdL1));
% plot(fC,  log10(psdC));
% plot(fR1, log10(psdR1));
% plot(fR2, log10(psdR2));
% legend("L2", "L1", "C", "R1", "R2")
% xlim([0 15]);


[H_l2l2, f_l2l2] = tfestimate(L2Flow2S(:,2), L2Flow2S(:,2), w, overlap, windowlen, fs); % calculate FRFs
[H_l2c, f_l2c]   = tfestimate(L2Flow2S(:,2), CFlow2S(:,2), w, overlap, windowlen, fs);
[H_l2r1, f_l2r1] = tfestimate(L2Flow2S(:,2), R1Flow2S(:,2), w, overlap, windowlen, fs);
[H_l2r2, f_l2r2] = tfestimate(L2Flow2S(:,2), R2Flow2S(:,2), w, overlap, windowlen, fs);
[H_l2l1, f_l2l1] = tfestimate(L2Flow2S(:,2), L1Flow2S(:,2), w, overlap, windowlen, fs);

i = 36;

mags = [abs(H_l2l2(i)), abs(H_l2l1(i)), abs(H_l2c(i)), abs(H_l2r1(i)), abs(H_l2r2(i))];
angles = [angle(H_l2l2(i)), angle(H_l2l1(i)), angle(H_l2c(i)), angle(H_l2r1(i)), angle(H_l2r2(i))];
xpos = [-226.31, -126.86, 0.00, 127.26, 226.58];
ypos = [-164.87, -249.61, -280.00, -249.41, -164.51];

t = 0:0.01:2*pi;
inputacc = sin(2*pi*3.5*t);
l2acc = mags(1)*sin(2*pi*3.5*t+angles(1));
l1acc = mags(2)*sin(2*pi*3.5*t+angles(2));
cacc = mags(3)*sin(2*pi*3.5*t+angles(3));
r1acc = mags(4)*sin(2*pi*3.5*t+angles(4));
r2acc = mags(5)*sin(2*pi*3.5*t+angles(5));

order = 4;
xvals = xpos(1):xpos(end);
for k = 1:length(t)
    accs = [l2acc(k), l1acc(k), cacc(k), r1acc(k), r2acc(k)];
    coefs = polyfit(xpos, accs, order);
    f = polyval(coefs, xvals);
    plot(xvals, f);
    hold on;
    plot(xvals, zeros(size(xvals)))
    ylim([-5 5])
    pause(0.1);
    clf
end




    

