%--------------------------------------------------------------------------
%	PROGRAM: Calculate the self-induction function by several methods
%                                                           
%	DEPENDENCY: 
%
%                                                         Zhangbo
%                                                        2021.11.14
%--------------------------------------------------------------------------
close all;clear all;clc

eps = 1e-6;
ka = (0:0.01:10);

% mode
m = 1;

% s = 1 for case in which the disturbance is retrograde,
% s = -1 for case in which the disturbance is co-grade
s = -1;

% Bessel function of first kind
myfun1 = @(x,mm) besselj(mm,x);

% the dispersion relation
myfun2 = @(ba,ka,s,m) 0.5*(besselj(m-1,ba)-besselj(m+1,ba))/besselj(m,ba)/ba- ...
                     0.5*(besselk(m-1,ka)+besselk(m+1,ka))/besselk(m,ka)/ka+ ...
                     s*m*sqrt(ba^2+ka^2)/ka/ba^2;

% find all roots of besseljm in the specified interval
x0 = [0 15]; % the specified interval
x_root = [];
r_old = -inf;
xstep = 1;
for i = x0(1):xstep:x0(2)
    r =fzero(@(x) myfun1(x,m),i);
    if abs(r-r_old) > 1e-5
        x_root = [x_root, r];
    end
    r_old = r;
end

% number of the considered root from the dispersion relation
N = length(x_root)-1;
omega = zeros(length(ka),N);

for j = 1:length(ka)
    for i = 1:N
        if ka(j) == 0
            % The solution at x = 0 is obtained by the limiting behavior
            % of the dispersion relation
            ba_root = x_root(i);
            
            if (i == 1) && (s ~= 1)
                % For the case s = -1, there is no solution in this interval
                omega(j,i) = 0;
            elseif i == 1
                omega(j,i) = 0;
            else
                omega(j,i) = 2*ka(j)*s/sqrt(ba_root^2+ka(j)^2)-m;
            end
        else
            xn0 = [x_root(i)+eps x_root(i+1)-eps];

            if (i == 1) && (s ~= 1)
                % For the case s = -1, there is no solution in this interval
                omega(j,i) = 0;
            else
                ba_root = fzero(@(ba) myfun2(ba,ka(j),s,m),xn0);
                
                omega(j,i) = 2*ka(j)*s/sqrt(ba_root^2+ka(j)^2)-m;
            end
        end
    end

end

% crow angular frequency: Gamma/(2*pi)*0.5*k^2*((cos(kd)-1)/kd^2+sin(kd)/kd-Ci(kd)) with d = 2*0.321*a
% that is: Gamma/(2*pi*a^2)*0.5*ka^2*((cos(kd)-1)/kd^2+sin(kd)/kd-Ci(kd)) with d = 2*0.321*a
% normalize: 0.5*ka^2*((cos(kd)-1)/kd^2+sin(kd)/kd-Ci(kd)) with d = 2*0.321*a
omega_crow = 0.5*ka.^2.*((cos(ka*0.642)-1)./(ka*0.642).^2+sin(ka*0.642)./(ka*0.642)-cosint(ka*0.642));
omega_LwAsymptotic = 0.5*ka.^2.*(log(2./ka)-0.577+1/4);
omega_Fit = ka.^2./(2+0.955*ka+0.438*ka.^2).*(log((2+2.151*ka)./ka)-0.577+1/4);

omega_crow = omega_crow';
omega_LwAsymptotic = omega_LwAsymptotic';
omega_Fit = omega_Fit';


if s==1
    p = plot(ka',omega,ka',omega_crow,'-.',ka',omega_LwAsymptotic,'--',ka',omega_Fit,':');
    %p(1).Marker = 'o';
    %p(1).MarkerIndices = 1:fix(length(ka)/25):length(ka);
    p(end).LineWidth = 1.5;
    axis([0 10 -1 1])
    grid on
    legend('exact:n=0','exact:n=1','exact:n=2','exact:n=3',...
           'crow:n=0','long-wavelength asymptotic:n=0','numerical fit:n=0',...
           'Location','Best')
else
    p = plot(ka',omega);
    axis([0 10 -3 1])
    grid on
    legend('exact:n=0','exact:n=1','exact:n=2','exact:n=3',...
           'Location','Best')
end
title('Vortex filament self-rotation rates for the bending modes with n radial node','Interpreter','latex')
xlabel('$ka$','Interpreter','latex')
ylabel('$\varpi_n/\Omega$','Interpreter','latex')
y0 = line([0,10],[0,0],'LineStyle','--','Color','red');
y0.Annotation.LegendInformation.IconDisplayStyle ='off';

fileID = fopen('Self-inductionFunction.dat','w');
if s==1
    for j=1:1:length(ka)
        fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\r\n',...
                ka(j),omega(j,1),omega(j,2),omega(j,3),omega(j,4),omega_crow(j,1),omega_LwAsymptotic(j,1),omega_Fit(j,1));
    end
else
    for j=1:1:length(ka)
        fprintf(fileID,'%6.2f %12.8f %12.8f %12.8f\r\n',...
                ka(j),omega(j,2),omega(j,3),omega(j,4));
    end
end