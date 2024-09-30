% author: Milos Nikolic, mnikolic@princeton.edu

clear;
addpath('auxiliary_matlab_functions')
load('rawProfiles_eve.mat')
%% extract data from data struct
Nemb=length(data);
L=[data.L];
age=[data.age];
thickness=[data.dist]; %cellularization membrane thickness

profiles_eve=nan(1000,Nemb);
y=(1:1000)/1000;
for k=1:Nemb
    profiles_eve(:,k)=data(k).Eve;
end
%% plot Fig.1 C
PROFILES=profiles_eve;
bw=5;
age_bins=bw/2:bw:60; %bin centers

mp=mean(PROFILES(:,80:end),2);
sp=nanstd(PROFILES(:,80:end),[],2);
mp_min=min(mp);
mp_max=max(mp);
cmap=colormap('turbo');
col_idx=round(age_bins./max(age_bins)*256);
kymograph_profiles=zeros(length(y),length(age_bins));

figure(1);
ax=axis();
set(gca,'fontname','Calibri','linewidth',1.2,'fontsize',14)
hold on
for k=1:length(age_bins)
    plot([0 1],age_bins(k)*[1 1],'-','linewidth',0.5,'color',0.8*[1 1 1])
end
for k=1:length(age_bins)
    curr_idx=find(age>age_bins(k)-bw/2 & age<age_bins(k)+bw/2);
    curr_profile=mean(PROFILES(:,curr_idx),2);
    curr_profile=(curr_profile-mp_min)./(mp_max-mp_min);
    plot(y,age_bins(k)+curr_profile*bw*1.2,'-','linewidth',2,'color',cmap(col_idx(k),:))
    kymograph_profiles(:,k)=curr_profile;
end

ylabel('Time in n.c. 14 (min)')
xlabel('x/L')
ylim([0 60+bw])
set(gca,'ytick',age_bins,'yticklabel',age_bins+bw/2)
title('Eve')
pbaspect([1,2,1])
hold off
%% remove embryos that have age t<30 min
data=data(71:end);
Nemb=length(data);
L=[data.L];
age=[data.age];
thickness=[data.dist];

profiles_eve=nan(1000,Nemb);
y=(1:1000)/1000;
for k=1:Nemb
    profiles_eve(:,k)=data(k).Eve;
end
%% prepare stripes for fitting
y=(1:1000)./1000;

%subtract a bakground value for the profiles for better fitting from a
%region where there is no expression.
bak=mean2(profiles_eve(50:250,:))+std2(profiles_eve(50:250,:)); 
profiles_eve_clean=profiles_eve-bak;
%

profiles_eve_clean(1:300,:)=nan;
profiles_eve_clean(850:1000,:)=nan;


%% inspect
figure(6)
plot(y,profiles_eve_clean)
xlabel('x_s')
ylabel('g')
%% fit eve
Nemb_clean=size(profiles_eve_clean,2);
pp_eve=nan(1,Nemb_clean); %peak position
all_pp_eve=nan(7,Nemb_clean);
all_param_eve=nan(21,Nemb_clean);
all_pp_f_eve=nan(7,Nemb_clean);

[pp,param,fval,pp_found]=fit_eve(y,profiles_eve_clean(:,1));
pp_eve(1)=pp;
all_param_eve(:,1)=param';
all_pp_eve(:,1)=(param(1:3:end)');
all_pp_f_eve(1:length(pp_found),1)=pp_found';

for k=2:Nemb_clean
    [pp,param,fval,pp_found]=fit_eve(y,profiles_eve_clean(:,k));
    pp_eve(k)=pp;
    all_param_eve(:,k)=param';
    all_pp_eve(:,k)=(param(1:3:end)');
    all_pp_f_eve(1:length(pp_found),k)=pp_found';
    title(num2str(k));
    pause(0.005);
    disp(['Fitted ',num2str(k), ' out of ',num2str(Nemb_clean),' expression profiles.'])
end

%% time correct eve

std_pp=nan(7,1);
mean_pp=nan(7,1);
std_pp_err=nan(7,1);
all_pp_eve_corr=nan(size(all_pp_eve));
t0=45; %set correction time to 45 min in nc14

%plotting options
def_colors=[0, 0.4470, 0.7410;...
0.8500, 0.3250, 0.0980;...
0.9290, 0.6940, 0.1250;...
0.4940, 0.1840, 0.5560;...
0.4660, 0.6740, 0.1880;...
0.3010, 0.7450, 0.9330;...
0.6350, 0.0780, 0.1840];
hf3=figure(3);
ax1=subplot(1,2,1);
hold(ax1,'on');
ax2=subplot(1,2,2);
hold(ax2,'on');

fit_age_coeff=nan(7,2);

for pp_idx=1:7
pp_curr_eve=all_pp_eve(pp_idx,:);

options = fitoptions('Method','NonlinearLeastSquares', 'Robust', 'LAR','StartPoint',[mean(pp_curr_eve),max(pp_curr_eve)-min(pp_curr_eve)]);
ft = fittype( 'poly1' );
[fitresult, gof] = fit( age', pp_curr_eve',ft, options);

pp_curr_eve_corr=pp_curr_eve-fitresult.p1*(age-t0);


fit_age_coeff(pp_idx,1)=fitresult.p1;
fit_age_coeff(pp_idx,1)=fitresult.p2;

all_pp_eve_corr(pp_idx,:)=pp_curr_eve_corr;

std_pp(pp_idx)=std(pp_curr_eve_corr);
mean_pp(pp_idx)=nanmean(pp_curr_eve_corr);
stats = bootstrp(1000,@(x)std(x),pp_curr_eve_corr);
std_pp_err(pp_idx)=std(stats);


plot(ax2,age,pp_curr_eve_corr,'.','color',def_colors(mod(pp_idx,8),:));
axis(ax2,[30 60 0.2 0.9]);

plot(ax1,age,pp_curr_eve,'.','color',def_colors(mod(pp_idx,8),:));
plot(ax1,[30 60],fitresult.p1*[30 60]+fitresult.p2,'-','color',def_colors(mod(pp_idx,8),:));
axis(ax1,[30 60 0.2 0.9]);

xlabel(ax1,'Embryo age (min)');
ylabel(ax1, 'Peak position x/L');
title(ax1,'Raw');

xlabel(ax2,'Embryo age (min)');
ylabel(ax2, 'Peak position x/L');
title(ax2,'Time corrected')
end

set(ax1,'fontsize',14);
set(ax2,'fontsize',14);

%% calculate eve precision

% calculate error based on anchoring. you only need the length variation.
x_fine=linspace(0,480,200);
sigma_A_fine=x_fine.*std(1./L);
sigma_P_fine=sqrt((mean(L)-x_fine).^2).*std(1./L);
sigma_AP_fine=sqrt(1./(1./sigma_A_fine.^2+1./sigma_P_fine.^2));
y_fine=x_fine./mean(L);

stats=bootstrp(1000,@(x) x_fine.*std(1./x),L);
% sigma_A_fine2=nanmean(stats);
sigma_A_fine_err=nanstd(stats);


stats=bootstrp(1000,@(x) sqrt((mean(L)-x_fine).^2).*std(1./x),L);
% sigma_P_fine2=nanmean(stats);
sigma_P_fine_err=nanstd(stats);


stats=bootstrp(1000,@(x) sqrt(1./(1./(x_fine.*std(1./x)).^2+1./(sqrt((mean(L)-x_fine).^2).*std(1./x)).^2)),L);
% sigma_AP_fine2=nanmean(stats);
sigma_AP_fine_err=nanstd(stats);
% plot precision sigma_x/L and anchoring
hf=figure(4);
f=1.5;
pbaspect([2 1 1])

def_colors=[0, 0.4470, 0.7410;...
0.8500, 0.3250, 0.0980;...
0.9290, 0.6940, 0.1250;...
0.4940, 0.1840, 0.5560;...
0.4660, 0.6740, 0.1880;...
0.3010, 0.7450, 0.9330;...
0.6350, 0.0780, 0.1840];
    

hold on

pgon1 = polyshape([y_fine fliplr(y_fine)],[sigma_A_fine+sigma_A_fine_err, fliplr(sigma_A_fine-sigma_A_fine_err)]);
hp_pgon1=plot(pgon1,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.5,'edgecolor','none');
hp4=plot(y_fine,sigma_A_fine,'k--');

pgon2 = polyshape([y_fine fliplr(y_fine)],[sigma_P_fine+sigma_P_fine_err, fliplr(sigma_P_fine-sigma_P_fine_err)]);
hp_pgon2=plot(pgon2,'FaceColor',[0.5 1 0.5],'FaceAlpha',0.5,'edgecolor','none');
hp5=plot(y_fine,sigma_P_fine,'k--');

pgon3 = polyshape([y_fine fliplr(y_fine)],[sigma_AP_fine+sigma_AP_fine_err, fliplr(sigma_AP_fine-sigma_AP_fine_err)]);
hp_pgon3=plot(pgon3,'FaceColor',[ 0.5 0.5 1 ],'FaceAlpha',0.25,'edgecolor','none');
hp6=plot(y_fine,sigma_AP_fine,'k-');
axis([ 0 1 0 0.02]);

load('Fig1_cephalicFurrow.mat'); %data from Liu et al. PNAS 2013. Collected data for 2XA reference fly line.
hp=errorbar(mean_pp,std_pp,std_pp_err,'b.');
errCF=nanstd(bootstrp(100,@(x) std(x),CF));
hpCF=errorbar(nanmean(CF),nanstd(CF),errCF*2,'kp','markersize',8,'markerfacecolor','k');
xlabel('x/L');
ylabel('\sigma_{x/L}');
legend([hp, hp_pgon1, hp_pgon2, hp_pgon3, hpCF],{'Eve','A anchor','P anchor','AP anchor','CF'},'location','northoutside');
set(gca,'fontsize',14,'fontname','helvetica','tickdir','in','linewidth',1);
grid on
box on
hold off
%% plot eve positions vs. L
figure(5)
hold on
plot([0 550],[0 550],'k--')
for k=1:7
    options = fitoptions('Method','NonlinearLeastSquares','Robust','LAR','StartPoint',[mean_pp(k)]);
    ft = fittype( 'a*x' );
    [fitresult, gof] = fit( L', all_pp_eve_corr(k,:)'.*L',ft);
    plot([0 550],fitresult.a*[0 550],'-','color',0.7*[1 1 1])
    pause(0.005);
end
plot(L,L'.*all_pp_eve_corr','.');


hold off
axis([0 550 0 550])
xlabel('x (\mum)')
ylabel('Absolute stripe position (\mum)');
title('Eve')