%% E101 HW6 MATLAB code
% Rafael Burger (rburger@hmc.edu) and Raja Batra (rbatra@hmc.edu)
% 11/29/2022

load("RVCradle.mat"); 

fs = 5000;     % sampling frequency of data in RVCradle
order = 3000;  % order for non-Fourier PSD estimation
N = 800000;    % section of data to analyze for non-Fourier PSD estimation

blocktime = 10;               % block time in seconds (to average FRF)
windowlen = blocktime * fs;   % convert block time from seconds to sample num
overlap = 0 * windowlen;      % block overlap for use in calculating FRFs
w = hann(windowlen);          % hanning window of specified windowlen for use calculating FRFs

%% Orientations and Positions of the Sensors
angles = [-144, -117, -90, -63, -36];
xcoords = [-226.31, -126.86, 0.00, 127.26, 226.58];
ycoords = [-164.87, -249.61, -280.00, -249.41, -164.51];

%% Calculate PSDs of all signals
% L2 = L2Flow2S(1:N,2);       % truncate data to allow faster PSD estimation
% L1 = L1Flow2S(1:N,2);
% C = CFlow2S(1:N,2);
% R1 = R1Flow2S(1:N,2);
% R2 = R2Flow2S(1:N,2);

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

%% Calculate FRFs of signals with L1 as reference
[H_l2, f_l2] = tfestimate(L1Flow2S, L2Flow2S, w, [], windowlen, fs); 
[H_c, f_c] = tfestimate(L1Flow2S, CFlow2S, w, [], windowlen, fs);
[H_r1, f_r1] = tfestimate(L1Flow2S, R1Flow2S, w, [], windowlen, fs);
[H_r2, f_r2] = tfestimate(L1Flow2S, R2Flow2S, w, [], windowlen, fs);

figure(1); 
hold on;
plot(f_l2, abs(H_l2(:,2)));
plot(f_c, abs(H_c(:,2)));
plot(f_r1, abs(H_r1(:,2)));
plot(f_r2, abs(H_r2(:,2)));
title("FRFs")
legend("L2", "C", "R1", "R2")
xlabel("frequency (hz)")
ylabel("Magnitude")
xlim([0 20])

%% Extract FRF magnitude and phase at 3.6Hz resonance from each FRF
ind   = 37;
magl2 = abs(H_l2(ind, 2));
magl1 = 1;  % mag = 1 bc l1 is the reference
magc  = abs(H_c(ind, 2));
magr1 = abs(H_r1(ind, 2));
magr2 = abs(H_r2(ind, 2));

phasel2 = angle(H_l2(ind, 2));
phasel1 = 0;  % phase = 0 bc l1 is the reference
phasec = angle(H_c(ind, 2));
phaser1 = angle(H_r1(ind, 2));
phaser2 = angle(H_r2(ind, 2));

%% Generate change in position vectors
t = 0:0.01:1/3.6; % time vector (10 oscillations)
w = 2*pi*3.6; % convert 3.6Hz frequency to radians 
% assuming arbitrary input acceleration signal at L1 = sin(w*t)
% Then output at each other point will be |H(jw)|sin(wt+angle(H(jw))
% double integrate this to get position = -(|H(jw)|/w^2)sin(wt+angle(H(jw))
% change in position vectors:
l2 = (-magl2/w^2)*sin(w*t+phasel2);
l1 = (-magl1/w^2)*sin(w*t+phasel1);
lc = (-magc/w^2)*sin(w*t+phasec);
r1 = (-magr1/w^2)*sin(w*t+phaser1);
r2 = (-magr2/w^2)*sin(w*t+phaser2);

%% Plot vibrating and original cradle for each timestep
figure(2);
while (1)
    for i = 1:length(t)
        m = [l2(i) l1(i) lc(i) r1(i) r2(i)]; %reformat dpos vectors
    
        % break dpos into x and y components based on angle of sensor
        dx = m .* cosd(angles); 
        dy = m .* sind(angles);
        
        amp = 500; %% amt by which to amplify the new signal
        
        %upate cradle from original using changes in x and y pos
        newxcoords = xcoords + amp*dx;
        newycoords = ycoords + amp*dy;
        
        % spline cradle positions onto a finer x-axis
        xq   = xcoords(1):0.1:xcoords(end);
        y    = interp1(xcoords, ycoords, xq, 'spline'); 
        newy = interp1(newxcoords, newycoords, xq, 'spline');    
        
        hold on;
        plot(xq, y, 'linewidth', 1);
        plot(xq, newy, 'linewidth', 1);
        ylim([-300 -160])
        pause(0.025)
        title("Cradle Position vs Time")
        subtitle(sprintf("Time: %d seconds", t(i)))
        ylabel("Position (?)")
        xlabel("Position (?)")
        legend("No input", "3.6Hz input")
        clf
    end
end


