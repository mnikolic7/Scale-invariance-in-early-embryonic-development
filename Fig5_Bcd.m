%This script plots the final results for calculation of delta I as described in
%SI Appendix, Sections B and D for Bicoid. It also reconstructs the Fig5 of
%the main text. 

%author: Milos Nikolic, mnikolic@princeton.edu


clear;
clc;
%% load Bcd profile data
% this data includes all profiles collected from Liu et al. PNAS 2013 for
% the reference fly line 2XA (see methods within). We excluded a small
% number of embryos (<30) that did not have complete measurement of Bcd
% concentration between 0.05 x/L to 0.90 x/L. 
% Intensity is normalized to the maximum of the mean expression profile
% (in the same way as for the gap genes).
load('Fig5_rawProfiles_Bcd.mat');
%% plot binned lengths (abs_coordinates) Fig5A
%
bw=10;
L_bins=440:bw:520; %bin centers

CP=profiles_bcd;
mp=nanmean(CP(:,1:end),2);
sp=nanstd(CP(:,1:end),[],2);
mp_min=min(mp);
mp_max=max(mp);
%set up colors for L bins
cmap=colormap('turbo');
col_idx=ceil(  ( L_bins-min(L_bins)+bw/2 )./(max(L_bins)-min(L_bins))*256  );
col_idx=[1 col_idx(1:end-1)];

binned_profiles=zeros(length(xs),length(L_bins)); 

hf=figure(1);
ax=subplot(1,2,1);
set(ax,'fontname','Helvetica','linewidth',1.2,'fontsize',12)
hold on

for k=2:length(L_bins)
    curr_idx=find(L>L_bins(k)-bw/2 & L<L_bins(k)+bw/2);
    curr_profile=nanmean(CP(:,curr_idx),2);
    curr_profile=(curr_profile-mp_min)./(mp_max-mp_min);
    curr_profile(82:end)=nan; %cut off profiles, since at posterior they are very close to zero and we are using log scale on the y axis.

    plot(ax,xs*L_bins(k),curr_profile,'-','linewidth',2,'color',cmap(col_idx(k),:))
    binned_profiles(:,k)=curr_profile;
end
ylabel('C_{bcd}(x)');
set(gca,'yscale','log')
xlabel(ax,'x (\mum)')

box on
xlim(ax,[0 500])
ylim([-10^-2 1]);
set(ax,'ytick',10.^[-2:1:0],'yticklabel',{'10^{-2}','10^{-1}','10^{0}'});
set(ax,'xtick',0:100:500)
hcb=colorbar;
ylabel(hcb,'L (\mum)','fontsize',14,'fontname','Helvetica'); 
pbaspect(ax,[1,1,1])
hold off
%% plot binned lenghts (rel_coordinates) Fig 5B
%
bw=10;
L_bins=440:bw:520; %bin centers
CP=profiles_bcd;
mp=nanmean(CP(:,1:end),2);
sp=nanstd(CP(:,1:end),[],2);
mp_min=min(mp);
mp_max=max(mp);
%set up colors for L bins
cmap=colormap('turbo');
col_idx=ceil(  ( L_bins-min(L_bins)+bw/2 )./(max(L_bins)-min(L_bins))*256  );
col_idx=[1 col_idx(1:end-1)];
binned_profiles=zeros(length(xs),length(L_bins)); 

hf=figure(1);
ax=subplot(1,2,2);
set(ax,'fontname','Helvetica','linewidth',1.2,'fontsize',12)
hold on
for k=2:length(L_bins)
    curr_idx=find(L>L_bins(k)-bw/2 & L<L_bins(k)+bw/2);
    curr_profile=nanmean(CP(:,curr_idx),2);
    curr_profile=(curr_profile-mp_min)./(mp_max-mp_min);
    curr_profile(82:end)=nan; %cut off profiles, since at posterior they are very close to zero and we are using log scale on the y axis.
    plot(ax,xs,curr_profile,'-','linewidth',2,'color',cmap(col_idx(k),:))
    binned_profiles(:,k)=curr_profile;

end

ylabel('C_{bcd}(x)');
set(gca,'yscale','log')
xlabel(ax,'x (\mum)')

box on
xlim(ax,[0 1])
ylim([-10^-2 1]);
set(ax,'ytick',10.^[-2:1:0],'yticklabel',{'10^{-2}','10^{-1}','10^{0}'});
hcb=colorbar;
ylabel(hcb,'L (\mum)','fontsize',14,'fontname','Helvetica'); 
colormap turbo
caxis([L_bins(1), L_bins(end)]);
pbaspect(ax,[1,1,1])
hold off

%% load dI calculation of Bcd
load('Fig5_Bcd.mat');

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

axis([-0.001 0.021 -0.025 0.5]) 
grid
set(gca,'fontsize',14,'fontname','Helvetica')
grid on
set(gca,'Box','off')
set(gca,'fontname','helvetica','fontsize',14,'tickdir','in')
