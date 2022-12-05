%% E101 Final Project Data Investigation
% Rafael Burger (rburger@hmc.edu) and Raja Batra (rbatra@hmc.edu)
% 12/4/2022
clear
load("E101_FP_DataFiles.mat");

time = Pier2nose_CS(:,1);
pier_acc = Pier2nose_CS(:,2);
spillway_acc = Spillway_atPier2_S(:,2);

fs = 1000;

%% Plot acceleration records
figure(1)
% plot(time, nose_acc); % plot acc vs time
plot(pier_acc); % plot acc vs sample num
xlabel("time (s)");
ylabel("acceleration (g)");
title("Acceleration at Nose");

figure(2)
% plot(time, spillway_acc); % plot acc vs time
plot(spillway_acc); % plot acc vs sample num
xlabel("time (s)");
ylabel("acceleration (g)");
title("Acceleration at Spillway");

%% Frequency Domain Analysis
order = 3000;

starts = [389974, 525030, 689249];
ends = [462835, 587833, 820458];

% preallocate cell arrays (b/c vectors will have diff lengths)
pier = cell(3, 1);
spillway = cell(3,1);
p_psd = cell(3,1);
p_f = cell(3,1);
s_psd = cell(3,1);
s_f = cell(3,1);

% select "on" segments and compute PSDs
for i = 1:3
    pier{i} = pier_acc(starts(i):ends(i));
    spillway{i} = spillway_acc(starts(i):ends(i));
    [p_psd{i}, p_f{i}] = pburg(pier{i}, order, length(pier{i}), fs);
    [s_psd{i}, s_f{i}] = pburg(spillway{i}, order, length(spillway{i}), fs);
end

% Plot PSDs
figure(3);
for i = 1:3
    subplot(3,1,i)
    semilogy(p_f{i}, p_psd{i});
    xlabel("frequency (Hz)");
    ylabel("power");
    title(sprintf("Pier: Spill %d", i));
    xlim([0 30]);
end

figure(4);
for i = 1:3
    subplot(3,1,i)
    semilogy(s_f{i}, s_psd{i});
    xlabel("frequency (Hz)");
    ylabel("power");
    title(sprintf("Spillway: Spill %d", i));
    xlim([0 30]);
end


