clc
clear all
close all

lambda   = linspace(450e-9, 1020e-9, 1000);
k0       = (2*pi)./lambda;
M        = 50;

L        = 20e-3;

a_mmf    = 25e-6;
a_NCF    = 62.5e-6;

Ncore_NCF = 1.4510;
Next_NCF  = 1;
ncore_MMF = 1.4512;
ncld_MMF  = 1.4316;

V_NCF = (2*pi./lambda) * a_NCF * sqrt(Ncore_NCF^2 - Next_NCF^2);
V_MMF = (2*pi./lambda) * a_mmf * sqrt(ncore_MMF^2 - ncld_MMF^2);

m_idx     = (1:M);
u_0m_vals = (m_idx - 0.25)*pi;
for it = 1:8
    u_0m_vals = u_0m_vals + besselj(0, u_0m_vals)./besselj(1, u_0m_vals);
end

w_s = a_mmf * ((2./V_MMF).^0.5 + 0.23*(V_MMF.^-1.5) + 18.01*(V_MMF.^-6));

TL_lambda = zeros(size(lambda));
for m = 1:M
    u_0m = u_0m_vals(m);
    numerator   = 2 .* (w_s./a_NCF).^2 .* exp( -((w_s./a_NCF).^2) .* (u_0m^2/2) );
    denominator = besselj(1, u_0m)^2;
    eta_0m = numerator ./ denominator;
    beta_0m = sqrt( k0.^2 * Ncore_NCF^2 - (u_0m^2 ./ a_NCF^2) );
    TL_lambda = TL_lambda + eta_0m .* exp(1i * beta_0m * L);
end

TL_lambda = 10 * log10( abs(TL_lambda).^2 );

figure;
plot(lambda*1e9, TL_lambda, 'LineWidth', 1.2);
xlim([820 1020])
xlabel('Wavelength [nm]', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('T_L(\lambda) [dB]', 'FontSize', 18, 'FontWeight', 'bold');
title(['Simulated Transmission for L = ', num2str(L*1e3), ' mm'], ...
      'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 3);
box on;
grid off;