function plotGACRCResult(inf,CRCPulse)
% ************************************************************************
%plotGACRCResult this function plots the best result inf for the 
% synthetic bipolar pulse CRCpulse
% ########################################################################
% INPUTS
% ########################################################################
% inf:      struct containing the best result
% CRCPulse: synthetic pulse
% ########################################################################
% Philipp Krumey
% 15-03-2021
% University of Duisburg-Essen
% *************************************************************************

%%calculate Fourier Coefficients and Fourier Series for synthetic bipolar
% pulse  
halfL=floor(length(CRCPulse)/2);
CRCPulseSinCoef = SinCoefs(CRCPulse,length(inf.bestDNA.sinCoef)-1,halfL);
CRCPulseCosCoef = CosCoefs(CRCPulse,length(inf.bestDNA.cosCoef)-1,halfL);
CRCPulsefSeries = fSeries(CRCPulseSinCoef,CRCPulseCosCoef,halfL);

%Plot Pulses
figure
plot(-halfL:halfL,CRCPulse,'LineWidth',2,'DisplayName','Original');
hold on
plot(-halfL:halfL,CRCPulsefSeries,'LineWidth',2,'DisplayName','Fourier Expansion');
plot(-halfL:halfL,inf.bestDNA.pulse,'LineWidth',2,'DisplayName','Retrieval');
hold off
title('bipolar Pulse with Retrieval')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xlabel('Length [nm]')
ylabel('Strain Amplitude')
legend

%Plot Sin Coefficients
figure
plot(CRCPulseSinCoef(2:end),'.','MarkerSize',24,'DisplayName','Original');
hold on
plot(inf.bestDNA.sinCoef(2:end),'.','MarkerSize',24,'DisplayName','Retrieval');
hold off
xlim([0.5 length(inf.bestDNA.sinCoef)-0.5])
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xlabel('Sine Coefficients')
ylabel('Strain Amplitude')
legend

%Plot Cos Coefficients
figure
plot(CRCPulseCosCoef,'.','MarkerSize',24,'DisplayName','Original');
hold on
plot(inf.bestDNA.cosCoef,'.','MarkerSize',24,'DisplayName','Retrieval');
hold off
xlim([0.5 length(inf.bestDNA.cosCoef)+0.5])
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xlabel('Cosine Coefficients')
ylabel('Strain Amplitude')
legend

%calculate rocking curves for best result and synthetic pulse
theta = inf.theta;                                                  %angles of rocking curve
time = [0 19 39 47 77];                                             %time points to be shown
CRCPulseCRC = norm2unp(calcCRC(CRCPulse,theta,time));               %calculate rocking curves for CRCPulse
CRCPulsefSeriesCRC = norm2unp(calcCRC(CRCPulsefSeries,theta,time)); %calculate rockings curves for Fourier Series of CRCPulse
retrievalCRC = norm2unp(calcCRC(inf.bestDNA.pulse,theta,time));     %calculate rocking corves for best result

%calculate fitness for best result and synthetic pulse
calcFitness(CRCPulsefSeriesCRC(:,5),CRCPulseCRC(:,5),theta,[-0.6 -0.1 0.1 0.6]) %calculate fitness of Fourier Series of CRCPulse for time[5]
calcFitness(retrievalCRC(:,5),CRCPulseCRC(:,5),theta,[-0.6 -0.1 0.1 0.6])       %calculate fitness of best result for time[5]

%plot rocking curves
figure
t = tiledlayout(2,2,'TileSpacing',"compact");

ax1 = nexttile;
plot(theta,CRCPulseCRC(:,2),'LineWidth',1.5,'DisplayName','Original')
hold on
plot(theta,CRCPulsefSeriesCRC(:,2),'LineWidth',1.5,'DisplayName','Fourier Expansion')
plot(theta,retrievalCRC(:,2),'LineWidth',1.5,'DisplayName','Retrieval')
hold off
title(append(num2str(time(2)),' ps'))
set(gca,'YScale','log')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xticklabels(ax1,{})
legend

ax2 = nexttile;
plot(theta,CRCPulseCRC(:,3),'LineWidth',1.5,'DisplayName','Original')
hold on
plot(theta,CRCPulsefSeriesCRC(:,3),'LineWidth',1.5,'DisplayName','Fourier Expansion')
plot(theta,retrievalCRC(:,3),'LineWidth',1.5,'DisplayName','Retrieval')
hold off
title(append(num2str(time(3)),' ps'))
set(gca,'YScale','log')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
xticklabels(ax2,{})
yticklabels(ax2,{})

ax3 = nexttile;
plot(theta,CRCPulseCRC(:,4),'LineWidth',1.5,'DisplayName','Original')
hold on
plot(theta,CRCPulsefSeriesCRC(:,4),'LineWidth',1.5,'DisplayName','Fourier Expansion')
plot(theta,retrievalCRC(:,4),'LineWidth',1.5,'DisplayName','Retrieval')
hold off
title(append(num2str(time(4)),' ps'))
set(gca,'YScale','log')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)

ax4 = nexttile;
plot(theta,CRCPulseCRC(:,5),'LineWidth',1.5,'DisplayName','Original')
hold on
plot(theta,CRCPulsefSeriesCRC(:,5),'LineWidth',1.5,'DisplayName','Fourier Expansion')
plot(theta,retrievalCRC(:,5),'LineWidth',1.5,'DisplayName','Retrieval')
hold off
title(append(num2str(time(5)),' ps'))
set(gca,'YScale','log')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)
yticklabels(ax4,{})

linkaxes([ax1,ax3],'x')
linkaxes([ax2,ax4],'x')
linkaxes([ax1,ax2],'y')
linkaxes([ax3,ax4],'y')
xlabel(t,'Theta [°]','fontweight','bold')
end

