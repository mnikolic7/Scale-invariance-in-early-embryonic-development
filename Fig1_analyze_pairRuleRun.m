% author: Milos Nikolic, mnikolic@princeton.edu

clear;
addpath('auxiliary_matlab_functions')
load('rawProfiles_run.mat')
%% extract data from data struct
Nemb=length(data);
L=[data.L];
age=[data.age];
thickness=[data.dist];

profiles_run=nan(1000,Nemb);
y=(1:1000)/1000;
for k=1:Nemb
    profiles_run(:,k)=data(k).Run;
end
%% prepare stripes for fitting
y=(1:1000)./1000;

%subtract a bakground value for the profiles for better fitting from a
%region where there is no expression.
bak=mean2(profiles_run(50:150,:))+std2(profiles_run(50:150,:)); 
profiles_run_clean=profiles_run-bak;
%

profiles_run_clean(1:300,:)=nan;
profiles_run_clean(900:1000,:)=nan;


%% inspect
figure(6)
plot(y,profiles_run_clean)
xlabel('x_s')
ylabel('g')
%% fit run
Nemb_clean=size(profiles_run_clean,2);
pp_run=nan(1,Nemb_clean); %peak position
all_pp_run=nan(7,Nemb_clean);
all_param_run=nan(21,Nemb_clean);
all_pp_f_run=nan(7,Nemb_clean);

[pp,param,fval,pp_found]=fit_eve(y,profiles_run_clean(:,1)); %use function fit_eve as runt is similar and requires no additional changes to the fitting parameters.
pp_run(1)=pp;
all_param_run(:,1)=param';
all_pp_run(:,1)=(param(1:3:end)');
all_pp_f_run(1:length(pp_found),1)=pp_found';

for k=2:Nemb_clean
    [pp,param,fval,pp_found]=fit_eve(y,profiles_run_clean(:,k));
    pp_run(k)=pp;
    all_param_run(:,k)=param';
    all_pp_run(:,k)=(param(1:3:end)');
    all_pp_f_run(1:length(pp_found),k)=pp_found';
    title(num2str(k));
    pause(0.005);% or use waitforbuttonpress;
    disp(['Fitted ',num2str(k), ' out of ',num2str(Nemb_clean),' expression profiles.'])
end
all_pp_run=sort(all_pp_run,1); %make sure the stripes are ordered after fitting.


%% time correct

std_pp=nan(7,1);
mean_pp=nan(7,1);
std_pp_err=nan(7,1);
all_pp_run_corr=nan(size(all_pp_run));
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
pp_curr_run=all_pp_run(pp_idx,:);

options = fitoptions('Method','NonlinearLeastSquares', 'Robust', 'LAR','StartPoint',[mean(pp_curr_run),max(pp_curr_run)-min(pp_curr_run)]);
ft = fittype( 'poly1' );
[fitresult, gof] = fit( age', pp_curr_run',ft, options);

pp_curr_run_corr=pp_curr_run-fitresult.p1*(age-t0);


fit_age_coeff(pp_idx,1)=fitresult.p1;
fit_age_coeff(pp_idx,1)=fitresult.p2;

all_pp_run_corr(pp_idx,:)=pp_curr_run_corr;

std_pp(pp_idx)=std(pp_curr_run_corr);
mean_pp(pp_idx)=nanmean(pp_curr_run_corr);
stats = bootstrp(1000,@(x)nanstd(x),pp_curr_run_corr);
std_pp_err(pp_idx)=std(stats);


plot(ax2,age,pp_curr_run_corr,'.','color',def_colors(mod(pp_idx,8),:));
axis(ax2,[30 60 0.2 0.9]);

plot(ax1,age,pp_curr_run,'.','color',def_colors(mod(pp_idx,8),:));
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

%% calculate precision

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

hp=errorbar(mean_pp,std_pp,std_pp_err,'b.');

xlabel('x/L');
ylabel('\sigma_{x/L}');
legend([hp, hp_pgon1, hp_pgon2, hp_pgon3],{'Run','A anchor','P anchor','AP anchor'},'location','northoutside');
set(gca,'fontsize',14,'fontname','helvetica','tickdir','in','linewidth',1);
grid on
box on
hold off
%% plot stripe positions vs. L
figure(5)
figure(5)
hold on
plot([0 550],[0 550],'k--')
for k=1:7
    options = fitoptions('Method','NonlinearLeastSquares','Robust','LAR','StartPoint',[mean_pp(k)]);
    ft = fittype( 'a*x' );
    [fitresult, gof] = fit( L', all_pp_run_corr(k,:)'.*L',ft);
    plot([0 550],fitresult.a*[0 550],'-','color',0.7*[1 1 1])
end
plot(L,L'.*all_pp_run_corr','.');
hold off
axis([0 550 0 550])
xlabel('x (\mu m)')
ylabel('Absolute stripe position (\mu m)');
title('Run')
