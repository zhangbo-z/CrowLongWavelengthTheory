%--------------------------------------------------------------------------
%	PROGRAM: Calculate the growth rate of long-wavelength instability by
%            Crow' theory
%                                                           
%	DEPENDENCY: 
%
%                                                         Zhangbo
%                                                        2021.11.12
%--------------------------------------------------------------------------
close all;clear all;clc

% dimensionless wavenumber, kb
beta = (0:0.01:6)';

% ratio of the equivalent core size and the separation distance, ae/b
ae_b = 0.3;

% dimensionless cutoff distance, kd
delta = 2*ae_b*0.321*beta;%0.2*beta;

% first mutual-induction function
chi = beta.*besselk(1,beta);

% second mutual-induction function
phi = beta.^2.*besselk(0,beta)+beta.*besselk(1,beta);

% self-induction function
omega = 0.5*((cos(delta)-1)./delta.^2+sin(delta)./delta-cosint(delta));

% dimensionless amplification rate
% symmetric mode
sigmaS = sqrt((1-phi+beta.^2.*omega).*(1+chi-beta.^2.*omega));
% antisymmetric mode
sigmaA = sqrt((1+phi+beta.^2.*omega).*(1-chi-beta.^2.*omega));

% extract the real part and imag part, respectively
selS_real = sigmaS == real(sigmaS);
sigmaS_real = sigmaS(selS_real);
selS_imag = sigmaS ~= real(sigmaS);
sigmaS_imga = sigmaS(selS_imag);

selA_real = sigmaA == real(sigmaA);
sigmaA_real = sigmaA(selA_real);
selA_imag = sigmaA ~= real(sigmaA);
sigmaA_imga = sigmaA(selA_imag);

% for plotting
sigmaS_real_plot = sigmaS;
sigmaS_real_plot(sigmaS_real_plot ~= real(sigmaS_real_plot)) = 0;

sigmaA_real_plot = sigmaA;
sigmaA_real_plot(sigmaA_real_plot ~= real(sigmaA_real_plot)) = 0;

% plot results
plot(beta,sigmaS_real_plot,'-',beta,sigmaA_real_plot,'--','LineWidth',1.5)
axis([0 6 0 1.5])
grid on
legend('Mode S','Mode A','Location','Best','Interpreter','latex')
title('Growth rate as a function of the normalized axial wavenumber,$\delta/\beta=0.1926$','Interpreter','latex')
xlabel('$\beta$','Interpreter','latex')
ylabel('$\sigma^*$','Interpreter','latex')
figure

plot(2*pi./beta,sigmaS_real_plot,'-',2*pi./beta,sigmaA_real_plot,'--','LineWidth',1.5)
axis([0 15 0 1.5])
grid on
legend('Mode S','Mode A','Location','Best','Interpreter','latex')
title('Growth rate as a function of the normalized axial wavelength, $a_e/b_0=0.3$','Interpreter','latex')
xlabel('$\lambda/b_0$','Interpreter','latex')
ylabel('$\sigma^*$','Interpreter','latex')

fileID1 = fopen('CrowTheory-Result-wavenumber.dat','w');
fileID2 = fopen('CrowTheory-Result-wavelength.dat','w');
for i = 1:length(beta)
    fprintf(fileID1,'%6.2f %12.8f %12.8f\n',beta(i),sigmaS_real_plot(i),sigmaA_real_plot(i));
    fprintf(fileID2,'%6.2f %12.8f %12.8f\n',2*pi./beta(i),sigmaS_real_plot(i),sigmaA_real_plot(i));
end
