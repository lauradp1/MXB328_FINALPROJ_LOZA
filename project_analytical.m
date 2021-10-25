close all
outflows = solution.outflows;
avg_outflow = abs(mean(outflows(:,:)));
t_vals = solution.t;
Kc = solution.constants.Kc;
%%
C0 = 1.78; %g/L;
D = 4e-9; %m^2/s
lambda = 3e-9; %s^-1
C_crit = 1e-6; %g\L
V = avg_outflow.*(1/(24*60*60)); %0.0108 - 0.0166
eta = linspace(0,65,50000);
eta = round(eta,3);
eta(1) = 0.001;
t = t_vals;
xc = 5;
%%



%contaminant = C0/2*(exp(f1*(1-gamma)*eta - x1^2)*erfcx(x1) + exp(f1*(1+gamma)*eta - x2^2)*erfcx(x2));

%plot(t,eta,contaminant)
contaminant = zeros(size(1:10:length(t),2), length(eta));
n = 1;
for i = 1:100:length(t)
    gamma  = (2*sqrt(D)*sqrt((V(i)^2./(4*D)) + lambda))/V(i);
    f2 = 2*sqrt(D);
    f1 = V(i)/(2*D);
    x1 = (eta - V(i)*gamma * t(i)*(24*60*60))/(f2*sqrt(t(i)*(24*60*60)));
    x2 = (eta + V(i)*gamma * t(i)*(24*60*60))/(f2*sqrt(t(i)*(24*60*60)));
    contaminant(n,:) = C0/2*(exp(f1*(1-gamma).*eta - x1.^2).*erfcx(x1) + exp(f1*(1+gamma).*eta - x2.^2).*erfcx(x2));
    figure(1)
    subplot(2,1,1)
    caption = sprintf('t = %0.2f (y)',(t(i)/(365)));
    plot(eta,contaminant(n,:),'LineWidth',3,'DisplayName',caption)
    %ylim([0 1.5e-6])
    hold on
    n = n+1;
end
legend show
figure(1)
xline(eta(38462),'-',{'Creek Location'},'FontSize',15,'LineWidth',3,'HandleVisibility','off')
%yline(1e-6,'-',{'Critical',' Concentration'},'FontSize',15,'LineWidth',3,'HandleVisibility','off')
title(['Contaminant Concentration Profile for $K_c =$ ',num2str(Kc),'$md^{-1}$'],'Interpreter', 'latex')
xlabel('Distance From Centre of Landfill $(m)$','Interpreter', 'latex')
ylabel('Contaminant Concentration $\left(gL^{-1}\right)$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
%%
contaminant = zeros(size(1:10:length(t),2), length(eta));
n = 1;
for i = 1:100:length(t)
    gamma  = (2*sqrt(D)*sqrt((V(i)^2./(4*D)) + lambda))/V(i);
    f2 = 2*sqrt(D);
    f1 = V(i)/(2*D);
    x1 = (eta - V(i)*gamma * t(i)*(24*60*60))/(f2*sqrt(t(i)*(24*60*60)));
    x2 = (eta + V(i)*gamma * t(i)*(24*60*60))/(f2*sqrt(t(i)*(24*60*60)));
    contaminant(n,:) = C0/2*(exp(f1*(1-gamma).*eta - x1.^2).*erfcx(x1) + exp(f1*(1+gamma).*eta - x2.^2).*erfcx(x2));
    figure(1)
    subplot(2,1,2)
    caption = sprintf('t = %0.2f (y)',(t(i)/(365)));
    plot(eta,contaminant(n,:),'LineWidth',3,'DisplayName',caption)
    ylim([0 1.5e-6])
    hold on
    n = n+1;
end
legend('Location', 'SouthEast')
figure(1)
xline(eta(38462),'-',{'Creek Location'},'FontSize',15,'LineWidth',3,'HandleVisibility','off')
yline(1e-6,'-',{'Critical',' Concentration'},'FontSize',15,'LineWidth',3,'HandleVisibility','off')
title(['Contaminant Concentration Profile for $K_c =$ ',num2str(Kc),'$md^{-1}$'],'Interpreter', 'latex')
xlabel('Distance From Centre of Landfill $(m)$','Interpreter', 'latex')
ylabel('Contaminant Concentration $\left(gL^{-1}\right)$', 'Interpreter', 'latex')
set(gca,'FontSize',15)

%legend.FontSize = 10;
%%

%$\hat{\psi}$','Interpreter','latex'
dt = t(2) - t(1);
z = find(eta == 50);
h = contaminant(:,z);
k = isnan(contaminant(:,z));
h(k) = C0;


%int = 1/xc * cumsum(contaminant(:,eta(z))*dt*V;

int = 1/xc * trapz(h*V);