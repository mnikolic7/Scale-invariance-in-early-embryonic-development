% author: Milos Nikolic, mnikolic@princeton.edu


clear;
addpath('auxiliary_matlab_functions')
load('rawProfiles_prd.mat')

%% extract data from data struct
Nemb=length(data);
L=[data.L];
age=[data.age];
thickness=[data.dist];

profiles_prd=nan(1000,Nemb);

y=(1:1000)/1000;
for k=1:Nemb
    profiles_prd(:,k)=data(k).Prd;
end
%% prepare stripes for fitting
y=(1:1000)./1000;

%subtract a bakground value for the profiles for better fitting from a
%region where there is no expression.
bak=mean2(profiles_prd(900:950,:))+std2(profiles_prd(900:950,:)); 
profiles_prd_clean=profiles_prd-bak;
%

profiles_prd_clean(1:300,:)=nan;
profiles_prd_clean(760:1000,:)=nan;


%% inspect
figure(6)
plot(y,profiles_prd_clean)
xlabel('x_s')
ylabel('g')

%% fit prd
% Fitting prd is slightly different becuase the stripe structure of prd
% is a bit more complicated, as the low expression domains between stripes
% tend to have some uneven expression which may appear as false stripes. 
% For more details please refer to:
% insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=1&ftext=FBgn0003145
% This does not prevent us from identifying the seven primary stripes of 
% prd, but the fitting procedure is adapted to account for this.
% 
% We manually provide the intial guesses for the peak positions here. 
% Note that along the dorsal profile used here, a small peak between 
% stripes 1 and 2 appears in the late stage of nc14. This is "peak #8" in 
% our fitting function. We ignore it in further analysis, as this is
% clearly an example of the uneven expression between stripes 1 and 2. 


Nemb_clean=size(profiles_prd_clean,2);

peak_position0=nan(8,Nemb_clean);
peak_position0(1,:)=0.34; %peak #1
peak_position0(2,1:end)=0.421; %peak #2
peak_position0(3,1:end)=0.477; %peak #3
peak_position0(4,1:end)=0.535; %peak #4
peak_position0(5,1:end)=0.594; %peak #5
peak_position0(6,1:end)=0.651; %peak #6
peak_position0(7,1:end)=0.709; %peak #7
peak_position0(8,1:end)=0.39; %peak #1.5 (Between 1 and 2 that appears late)

pp_prd=ones(size(peak_position0))*-1;
pp_prd(isnan(peak_position0))=nan;
all_pp_prd=pp_prd;
all_param_prd=-1*ones(30,Nemb_clean);

idx_curr=find(~isnan(peak_position0(:,1)));
pp0_curr=peak_position0(idx_curr,1);
[pp,param,fval]=fit_prd(y,profiles_prd_clean(:,1),pp0_curr);
pp_prd(idx_curr,1)=pp;
all_param_prd(3*idx_curr-2,1)=param(1:3:end)';
all_param_prd(3*idx_curr-1,1)=param(2:3:end)';
all_param_prd(3*idx_curr,1)=param(3:3:end)';
all_pp_prd(idx_curr,1)=(param(1:3:end)');

for k=2:Nemb_clean
    idx_curr=find(~isnan(peak_position0(:,k)));
    pp0_curr=peak_position0(idx_curr,k);
    [pp,param,fval]=fit_prd(y,profiles_prd_clean(:,k),pp0_curr);
    pp_prd(idx_curr,k)=pp;
    all_param_prd(3*idx_curr-2,k)=param(1:3:end)';
    all_param_prd(3*idx_curr-1,k)=param(2:3:end)';
    all_param_prd(3*idx_curr,k)=param(3:3:end)';
    all_pp_prd(idx_curr,k)=(param(1:3:end)');

    disp(k);
    pause(0.005);

end
%
all_pp_prd=sort(all_pp_prd,1);
all_pp_prd(2,:)=[]; 
%% time correct

std_pp=nan(7,1);
mean_pp=nan(7,1);
std_pp_err=nan(7,1);
all_pp_prd_corr=nan(size(all_pp_prd));
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
pp_curr_prd=all_pp_prd(pp_idx,:);

options = fitoptions('Method','NonlinearLeastSquares', 'Robust', 'LAR','StartPoint',[mean(pp_curr_prd),max(pp_curr_prd)-min(pp_curr_prd)]);
ft = fittype( 'poly1' );
[fitresult, gof] = fit( age', pp_curr_prd',ft, options);

pp_curr_prd_corr=pp_curr_prd-fitresult.p1*(age-t0);


fit_age_coeff(pp_idx,1)=fitresult.p1;
fit_age_coeff(pp_idx,1)=fitresult.p2;

all_pp_prd_corr(pp_idx,:)=pp_curr_prd_corr;

std_pp(pp_idx)=std(pp_curr_prd_corr);
mean_pp(pp_idx)=nanmean(pp_curr_prd_corr);
stats = bootstrp(1000,@(x)nanstd(x),pp_curr_prd_corr);
std_pp_err(pp_idx)=std(stats);


plot(ax2,age,pp_curr_prd_corr,'.','color',def_colors(mod(pp_idx,8),:));
axis(ax2,[30 60 0.2 0.9]);

plot(ax1,age,pp_curr_prd,'.','color',def_colors(mod(pp_idx,8),:));
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
legend([hp, hp_pgon1, hp_pgon2, hp_pgon3],{'Prd','A anchor','P anchor','AP anchor'},'location','northoutside');
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
    [fitresult, gof] = fit( L', all_pp_prd_corr(k,:)'.*L',ft);
    plot([0 550],fitresult.a*[0 550],'-','color',0.7*[1 1 1])
end
plot(L,L'.*all_pp_prd_corr','.');

hold off
axis([0 550 0 550])
xlabel('x (\mu m)')
ylabel('Absolute stripe position (\mu m)');
