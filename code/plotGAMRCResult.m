function plotGAMRCResult(inf,ModelPulse,MRC)

halfL=floor(length(ModelPulse)/2);
ModelPulseSinCoef = SinCoefs(ModelPulse,length(inf.bestDNA.sinCoef)-1,halfL);
ModelPulseCosCoef = CosCoefs(ModelPulse,length(inf.bestDNA.cosCoef)-1,halfL);
ModelPulsefSeries = fSeries(ModelPulseSinCoef,ModelPulseCosCoef,halfL);

%Plot Pulses
figure
plot(-halfL:halfL,ModelPulse,'LineWidth',2,'DisplayName','Model');
hold on
plot(-halfL:halfL,ModelPulsefSeries,'LineWidth',2,'DisplayName','Fourier Expansion');
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
plot(ModelPulseSinCoef(2:end),'.','MarkerSize',24,'DisplayName','Model');
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
plot(ModelPulseCosCoef,'.','MarkerSize',24,'DisplayName','Model');
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

%Plot Rocking Curves
theta = inf.theta;
time = [0 19 39 47 77];
ModelPulseCRC = norm2unp(calcCRC(ModelPulse,theta,time));
retrievalCRC = norm2unp(calcCRC(inf.bestDNA.pulse,theta,time));

calcFitness(ModelPulseCRC(:,5),MRC(:,5),theta,[-0.6 -0.1 0.1 0.6])
calcFitness(retrievalCRC(:,5),MRC(:,5),theta,[-0.6 -0.1 0.1 0.6])

figure
t = tiledlayout(2,2,'TileSpacing',"compact");

ax1 = nexttile;
plot(theta,ModelPulseCRC(:,2),'LineWidth',1.5,'DisplayName','Model')
hold on
plot(theta,MRC(:,2),'LineWidth',1.5,'DisplayName','Original')
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
plot(theta,ModelPulseCRC(:,3),'LineWidth',1.5,'DisplayName','Model')
hold on
plot(theta,MRC(:,3),'LineWidth',1.5,'DisplayName','Orignal')
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
plot(theta,ModelPulseCRC(:,4),'LineWidth',1.5,'DisplayName','Model')
hold on
plot(theta,MRC(:,4),'LineWidth',1.5,'DisplayName','Original')
plot(theta,retrievalCRC(:,4),'LineWidth',1.5,'DisplayName','Retrieval')
hold off
title(append(num2str(time(4)),' ps'))
set(gca,'YScale','log')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
set(gca,'linewidth',2)

ax4 = nexttile;
plot(theta,ModelPulseCRC(:,5),'LineWidth',1.5,'DisplayName','Model')
hold on
plot(theta,MRC(:,5),'LineWidth',1.5,'DisplayName','Original')
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

