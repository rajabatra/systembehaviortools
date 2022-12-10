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
order = 650;

[p_psd, p_f] = pburg(pier_acc, order, length(pier_acc), fs);
[sw_psd, sw_f] = pburg(spillway_acc, order, length(spillway_acc), fs);

% 
% 
% % Plot PSDs
% figure(3);
% semilogy(p_f, p_psd);
% xlabel("frequency (Hz)");
% ylabel("power");
% title("Pier PSD");
% xlim([0 30]);
% 
% 
% figure(4);
% semilogy(sw_f, sw_psd);
% xlabel("frequency (Hz)");
% ylabel("power");
% title("Spillway PSD");
% xlim([0 30]);


wn_p = [5.234, 7.37];
z_p = 0.0782;

cpmp = 2*z_p*wn_p; %%2 component vectors
kpmp = wn_p.^2;

wn_sw = [1.93, 3.276, 4.778, 5.957, 7.343];
z_sw = 0.045;

cswmsw = 2*z_sw*wn_sw; %% 4 component vectors
kswmsw = wn_sw.^2;

r = 10^2;
step = 0.01;
w = step:step:50;
x = zeros(size(w));

for p = 1:2
    for sw = 1:4
    a = -(w.^2) + 1i.*w.*cpmp(p) + kpmp(p);
    b = -1i.*w*cpmp(p)-kpmp(p);
    c = b;
    d = -(w.^2)*r+1i.*2.*cpmp(p) + 1i.*w.*cswmsw(sw)*r + kpmp(p) + kswmsw(sw)*r;
    f = (1i.*w).^-1;
    x = x + (a.*f)./(a.*d-b.*c);
    end
end

model_psd = -(w.^2).*x;

figure(3)
hold on;
yyaxis left
semilogy(w, abs(model_psd))
yyaxis right
semilogy(sw_f, sw_psd);
xlim([0 10]);
xlabel("frequency (hz)");
title("Real and Predicted PSD");
legend("Model PSD", "Real PSD");



