function [f_F,f_NF, PSD_F, PSD_NF] = e101hw5(x, fs, order)
N = length(x);
[Xjw, f_F] = fdomain(x, fs);
PSD_F = (N/fs)*(Xjw.*(conj(Xjw)));
[PSD_NF, f_NF] = pburg(x, order, N, fs);
%PSD_NF = 0;

% rms_F = rms(PSD_F);
% rms_NF = rms(PSD_NF);
% rms_error = 100*abs(rms_F-rms_NF)/rms_F;

% figure(1)
% semilogy(f_F,PSD_F);
% hold on;
% semilogy(f_NF,PSD_NF, 'LineWidth', 1.5);
% xlabel("frequency");
% ylabel("power");
% title("PSD: Fourier and Non-Fourier Method")
% % subtitle(sprintf("Pburg RMS Percent Error: %.1f %", rms_error))
% legend("fourier", sprintf("pburg (Order: %d)", order))
end
