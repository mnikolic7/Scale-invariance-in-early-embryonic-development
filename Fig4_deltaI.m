%This script plots the final results for calculation of delta I as described in
%SI Appendix, Sections B and D. 

%author: Milos Nikolic, mnikolic@princeton.edu
clear;
clc;
%% load dI calculation of Hb
%note that we calculate dI from 100 random subsamples of the whole dataset.
%variable NN indicates the sizes of the subsamples.
load('ready_for_submission\Fig4_Hb.mat');

% plot result
figure(1)
hold on
for i=1:100
    hplt_pts=scatter(1./NN,dI(i,:), 'Marker', 'o', 'markeredgecolor','none','markerfacecolor','k', 'MarkerfaceAlpha', 0.5);
    for k=1:length(hplt_pts)
        hplt_pts(k).SizeData = hplt_pts(k).SizeData./9;
    end
end

hold on
errorbar(1./NN,mean(dI),std(dI),'o','color',def_colors(1,:),'markerfacecolor',def_colors(1,:))
plot((0:0.001:0.05),AAdI2(:,1)*(0:0.001:0.05).^2 + AAdI2(:,2)*(0:0.001:0.05) + AAdI2(:,3),'k-','linewidth', 1);
errorbar(0,dI_at_inf,err,'o','color',def_colors(7,:),'markerfacecolor','none')
hold off
xlabel('$1/N_{\rm{em}}$','interpreter','latex')
ylabel('$\Delta I (\rm{bits})$','interpreter','latex')

axis([-0.001 0.021 -0.025 0.3]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
%% load dI calculation of Gt
load('ready_for_submission\Fig4_Gt.mat');

% plot result
figure(2)
hold on
for i=1:100
    hplt_pts=scatter(1./NN,dI(i,:), 'Marker', 'o', 'markeredgecolor','none','markerfacecolor','k', 'MarkerfaceAlpha', 0.5);
    for k=1:length(hplt_pts)
        hplt_pts(k).SizeData = hplt_pts(k).SizeData./9;
    end
end

hold on
errorbar(1./NN,mean(dI),std(dI),'o','color',def_colors(1,:),'markerfacecolor',def_colors(1,:))
plot((0:0.001:0.05),AAdI2(:,1)*(0:0.001:0.05).^2 + AAdI2(:,2)*(0:0.001:0.05) + AAdI2(:,3),'k-','linewidth', 1);
errorbar(0,dI_at_inf,err,'o','color',def_colors(7,:),'markerfacecolor','none')
hold off
xlabel('$1/N_{\rm{em}}$','interpreter','latex')
ylabel('$\Delta I (\rm{bits})$','interpreter','latex')

axis([-0.001 0.021 -0.025 0.3]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
%% load dI calculation of Kni

%In the middle part of Kni from 0.15 to 0.45 the small middle/anterior 
%stripe of Kni is expressed. This stripe varies systematically along the 
%dorsoventral axis. Note that these positions are excluded from this 
%analysis, as we cannot exclude that the slight error in DV orientation of 
%the embryos during imaging does not correlate with length,
%and would thus create artifacts.
% See:
% https://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=1&ftext=FBgn0001320
% for an example of the Kni expression profile.
load('ready_for_submission\Fig4_Hb.mat');

% plot result
figure(3)
hold on
for i=1:100
    hplt_pts=scatter(1./NN,dI(i,:), 'Marker', 'o', 'markeredgecolor','none','markerfacecolor','k', 'MarkerfaceAlpha', 0.5);
    for k=1:length(hplt_pts)
        hplt_pts(k).SizeData = hplt_pts(k).SizeData./9;
    end
end

hold on
errorbar(1./NN,mean(dI),std(dI),'o','color',def_colors(1,:),'markerfacecolor',def_colors(1,:))
plot((0:0.001:0.05),AAdI2(:,1)*(0:0.001:0.05).^2 + AAdI2(:,2)*(0:0.001:0.05) + AAdI2(:,3),'k-','linewidth', 1);
errorbar(0,dI_at_inf,err,'o','color',def_colors(7,:),'markerfacecolor','none')
hold off
xlabel('$1/N_{\rm{em}}$','interpreter','latex')
ylabel('$\Delta I (\rm{bits})$','interpreter','latex')

axis([-0.001 0.021 -0.025 0.3]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
%% load dI calculation of Kr
load('ready_for_submission\Fig4_Hb.mat');

% plot result
figure(4)
hold on
for i=1:100
    hplt_pts=scatter(1./NN,dI(i,:), 'Marker', 'o', 'markeredgecolor','none','markerfacecolor','k', 'MarkerfaceAlpha', 0.5);
    for k=1:length(hplt_pts)
        hplt_pts(k).SizeData = hplt_pts(k).SizeData./9;
    end
end

hold on
errorbar(1./NN,mean(dI),std(dI),'o','color',def_colors(1,:),'markerfacecolor',def_colors(1,:))
plot((0:0.001:0.05),AAdI2(:,1)*(0:0.001:0.05).^2 + AAdI2(:,2)*(0:0.001:0.05) + AAdI2(:,3),'k-','linewidth', 1);
errorbar(0,dI_at_inf,err,'o','color',def_colors(7,:),'markerfacecolor','none')
hold off
xlabel('$1/N_{\rm{em}}$','interpreter','latex')
ylabel('$\Delta I (\rm{bits})$','interpreter','latex')

axis([-0.001 0.021 -0.025 0.3]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
%% load dI calculation of combined gap genes
load('ready_for_submission\Fig4_gapCombined_HbGtKniKr.mat');

% plot result
figure(5)
hold on
for i=1:100
    hplt_pts=scatter(1./NN,dI(i,:), 'Marker', 'o', 'markeredgecolor','none','markerfacecolor','k', 'MarkerfaceAlpha', 0.5);
    for k=1:length(hplt_pts)
        hplt_pts(k).SizeData = hplt_pts(k).SizeData./9;
    end
end

hold on
errorbar(1./NN,mean(dI),std(dI),'o','color',def_colors(1,:),'markerfacecolor',def_colors(1,:))
plot((0:0.001:0.05),AAdI2(:,1)*(0:0.001:0.05).^2 + AAdI2(:,2)*(0:0.001:0.05) + AAdI2(:,3),'k-','linewidth', 1);
errorbar(0,dI_at_inf,err,'o','color',def_colors(7,:),'markerfacecolor','none')
hold off
xlabel('$1/N_{\rm{em}}$','interpreter','latex')
ylabel('$\Delta I (\rm{bits})$','interpreter','latex')

axis([-0.001 0.011 -0.025 1]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
