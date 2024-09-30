% load gap gene data and fit domain boundaries
% domain boundary is the midpoint between high and low intensity domains of
% expression
%author: Milos Nikolic, mnikolic@princeton.edu

clear;
clc;
addpath('auxiliary_matlab_functions')
%%
load('rawProfiles_gapGenes_Hb_Gt_Kni_Kr.mat') 
%gene expression is already normalized to the range of min and max of the
%mean expression profile within individual experimental sessions (see SI,
%Appendix, section D).
%%
Nemb=length(data);

age=[data.age];
L=[data.L];
profiles_hb=horzcat(data.Hb);
profiles_gt=horzcat(data.Gt);
profiles_kni=horzcat(data.Kni);
profiles_kr=horzcat(data.Kr);
xs=data(1).xs;
%% Fit all boundaries. 
%% hunchback boundary at ~0.47
disp 'fitting hb at 0.47'
Nemb=size(profiles_hb,2);
pp_hb1=nan(1,Nemb);
all_param_hb1=nan(4,Nemb);


for k=1:Nemb
    curr_profile=profiles_hb(:,k);
    curr_profile=curr_profile(200:630);
    xs_curr=xs(200:630);
    [pp,param]=fit_boundaryR(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(curr_profile)],'k--');
    axis([0 1 -0.1 2])
    grid on
    hold off
    pp_hb1(k)=pp;
    all_param_hb1(:,k)=param';
    title(num2str(k));
    pause(0.01);
end

%% hunchback posterior boundaries peak near 0.8 
%posterior peak of hb expression near 0.8 is weak very early in nc14.
%some very early embryos have very low expression, which makes it difficult
%to identify the boundaries. Ignore those data points. 
disp 'fitting hb at 0.8'
Nemb=size(profiles_hb,2);
pp_hb2=nan(2,Nemb);
all_param_hb2=nan(8,Nemb);
min_max_ratio=zeros(1,Nemb);
ratio_first_to_peak=zeros(Nemb,1);
for k=1:Nemb
    curr_profile=profiles_hb(:,k);
    curr_profile=curr_profile(650:900);
    xs_curr=xs(650:900);
    ratio_first_to_peak(k)=(max(curr_profile)-curr_profile(1)); %first point is the background in this case
    if ratio_first_to_peak(k)<0.05
        continue;
    end
    [pp,param]=fit_peak_splitLRmidpoint(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(curr_profile)],'k--');
    plot(param(5)*[1 1],[0, max(curr_profile)],'k--');
    hold off
    pp_hb2(:,k)=pp;
    all_param_hb2(:,k)=param';
    title(num2str(k));
    pause(0.01);
end

%% collect hb fits
pp_hb=[pp_hb1;pp_hb2];
all_param_hb=[all_param_hb1;all_param_hb2];

%% gt boundary at ~0.4 
disp 'fitting gt at 0.4'
Nemb=size(profiles_gt,2);
pp_gt1=nan(1,Nemb);
all_param_gt1=nan(4,Nemb);


for k=1:Nemb
    curr_profile=profiles_gt(:,k);
    curr_profile=curr_profile(250:520);
    xs_curr=xs(250:520);
    min_max_ratio(k)=max(curr_profile)./min(curr_profile);
    [pp,param]=fit_boundaryR(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1], [0 max(curr_profile)],'k--')
    axis([0 1 -0.1 1.5])
    grid on
    hold off
    pp_gt1(k)=pp;
    all_param_gt1(:,k)=param';
    title(num2str(k));
    pause(0.01);
end
%% gt boundary at ~0.22
disp 'fitting gt at 0.22'
Nemb=size(profiles_gt,2);
pp_gt2=nan(1,Nemb);
all_param_gt2=nan(4,Nemb);

for k=1:Nemb
    curr_profile=profiles_gt(:,k);
    curr_profile=curr_profile(70:500);
    xs_curr=xs(70:500);
    min_max_ratio(k)=max(curr_profile)./min(curr_profile);
    [pp,param]=fit_boundaryL(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1], [0 max(curr_profile)],'k--')
    axis([0 1 -0.1 1.5])
    grid on
    hold off

    pp_gt2(k)=pp;
    all_param_gt2(:,k)=param';
    title(num2str(k));
    pause(0.01);
end
%% gt peak at posterior
disp 'fitting gt at posterior'
Nemb=size(profiles_gt,2);
pp_gt3=nan(2,Nemb);
all_param_gt3=nan(8,Nemb);

for k=1:Nemb
    curr_profile=profiles_gt(:,k);
    curr_profile=curr_profile(550:900);
    xs_curr=xs(550:900);
    [pp,param]=fit_peak_splitLRmidpoint(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(curr_profile)],'k--');
    plot(param(5)*[1 1],[0, max(curr_profile)],'k--');
    axis([0 1 -0.1 1.5])
    grid on
    hold off
    pp_gt3(:,k)=pp;
    all_param_gt3(:,k)=param';
    title(num2str(k));
    pause(0.01);
end
%% collect gt boundaries
pp_gt=[pp_gt1;pp_gt2;pp_gt3];
all_param_gt=[all_param_gt1;all_param_gt2;all_param_gt2];

%% kni
disp 'fitting kni peak'
Nemb=size(profiles_kni,2);
pp_kni1=nan(2,Nemb);
all_param_kni1=nan(8,Nemb);

for k=1:Nemb
    curr_profile=profiles_kni(:,k);
    curr_profile=curr_profile(350:900);
    xs_curr=xs(350:900);
    [pp,param]=fit_peak_splitLRmidpoint(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(curr_profile)],'k--');
    plot(param(5)*[1 1],[0, max(curr_profile)],'k--');
    axis([0 1 -0.1 1.5])
    grid on
    hold off
    pp_kni1(:,k)=pp;
    all_param_kni1(:,k)=param';

    title(num2str(k));
    pause(0.0001); 
%     waitforbuttonpress;
end
%% knirps small peak
% this small peak is expressed only late in nc14.
% Fit only profiles that have this profile clearly visible. 
disp 'fitting kni small tiny peak'
Nemb=size(profiles_kni,2);
pp_kni2=nan(1,Nemb);
all_param_kni2=nan(4,Nemb);
pp_kni_x2=nan(1,Nemb);
all_param_kni_x2=nan(4,Nemb);
min_max_ratio=zeros(1,Nemb);
%
for k=1:Nemb
    curr_profile=profiles_kni(:,k);
    curr_profile=curr_profile(200:400);
    xs_curr=xs(200:400);
    min_max_ratio(k)=max(curr_profile)./min(curr_profile);
    if min_max_ratio(k)<3 && age(k)<60
        continue;
    end
    [pp,param]=fit_peak_simple(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(profiles_kni(:,k))],'k--');
    plot(xs_curr,param(3)/(param(2)*sqrt(2*pi))*exp(-0.5*(xs_curr-param(1)).^2./(param(2).^2))+param(4),'r-')
    axis([0 1 -0.1 1.5])
    grid on
    hold off
    pp_kni2(k)=pp;
    all_param_kni2(:,k)=param';
    title(num2str(k));
    pause(0.01); 

end
%%
disp 'fitting kni anterior position'
Nemb=size(profiles_kni,2);
pp_kni3=nan(1,Nemb);
all_param_kni3=nan(4,Nemb);

for k=1:Nemb
    curr_profile=profiles_kni(:,k);
    curr_profile=curr_profile(1:200);
    xs_curr=xs(1:200);
    [pp,param]=fit_boundaryR(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1], [0 max(curr_profile)],'k--')
    axis([0 1 -0.1 1.5])
    grid on
    hold off

 
    pp_kni3(k)=pp;
    all_param_kni3(:,k)=param';

    title(num2str(k));
    pause(0.001);
end
%% aggregate knirps results
pp_kni=[pp_kni1;pp_kni2;pp_kni3];
all_param_kni=[all_param_kni1;all_param_kni2;all_param_kni2];

%% fit kr
disp 'fitting kr, almost done.'
Nemb=size(profiles_kr,2);
pp_kr=nan(2,Nemb);
all_param_kr=nan(8,Nemb);

for k=1:Nemb
    curr_profile=profiles_kr(:,k);
    curr_profile=curr_profile(150:700);
    xs_curr=xs(150:700);
    min_max_ratio(k)=abs((max(curr_profile)-min(curr_profile))./min(curr_profile));
    if min_max_ratio(k)<5
        continue;
    end
    [pp,param]=fit_peak_splitLRmidpoint(xs_curr,curr_profile);
    plot(xs_curr,curr_profile,'o')
    hold on
    plot(param(1)*[1 1],[0, max(curr_profile)],'k--');
    plot(param(5)*[1 1],[0, max(curr_profile)],'k--');
    axis([0 1 -0.1 1.5])
    grid on
    hold off
    pp_kr(:,k)=pp;
    all_param_kr(:,k)=param;
    title(num2str(k));
    pause(0.01); 
end
%% plot kr for inspection
figure(1);
imagesc(1:Nemb,xs,(profiles_kr))
hold on
plot(1:Nemb,pp_kr','r-')
hold off
xlabel('Embryo #');
ylabel('x_s');

%%
% age=age_clean;
% L=L_clean;
if isrow(age)
    age=age';
end
if isrow(L)
    L=L';
end
%% a small number of profiles are outliers. Typically this is the case for 
% the very early embryos (<20min) when not all expression domains are
% defined well. This is most exaggerated for the 
% posterior boundary of giant, for example. 
% Those are either experimental errors or failures of the fitting
% procedure. 
% An effective way to remove them is the rmoutliers function
% from matlab which removes points >3 sigma away from the mean*. Otherwise,
% they can be manually identified, and provide the same result. If they are
% not accounted for, they greatly skew the estimate of the standard
% deviation.

%*rmoutliers uses average median deviation (MAD), which is more robust, but
%the function is scaled in such a way to remove the otliers beyond 3sigma
%(99.7% confidence interval).
%% correct for time, and remove outliers.
hb=pp_hb';
t0=45; %
Nb=size(hb,2);
hb_corr=nan(size(hb));

IDX_REMOVE=cell(4,4); %keep track of which points were outliers
for k=1:Nb
    idx_keep=~isnan(hb(:,k));
    fit_param=polyfit(age(idx_keep),hb(idx_keep,k),1);
    hb_corr(:,k)=hb(:,k)-(fit_param(1)*(age-t0));
    [~,idx_remove]=rmoutliers(hb_corr(:,k),'mean');
    hb_corr(idx_remove,k)=nan;
    IDX_REMOVE{1,k}=idx_remove;
end
hb_corr_x=hb_corr.*L;
% plot(age,hb_corr,'.'); %inspect
%%
gt=pp_gt';
Nb=size(gt,2);
gt_corr=nan(size(gt));
for k=1:Nb
    idx_keep=~isnan(gt(:,k));
    fit_param=polyfit(age(idx_keep),gt(idx_keep,k),1);
    gt_corr(:,k)=gt(:,k)-(fit_param(1)*(age-t0));
    [~,idx_remove]=rmoutliers(gt_corr(:,k),'mean');
    gt_corr(idx_remove,k)=nan;
    IDX_REMOVE{2,k}=idx_remove;
end
gt_corr_x=gt_corr.*L;
% plot(age,gt_corr,'.'); %inspect

kni=pp_kni';
Nb=size(kni,2);
kni_corr=nan(size(kni));
for k=1:Nb
    idx_keep=~isnan(kni(:,k));
    fit_param=polyfit(age(idx_keep),kni(idx_keep,k),1);
    kni_corr(:,k)=kni(:,k)-(fit_param(1)*(age-t0));
    [~,idx_remove]=rmoutliers(kni_corr(:,k),'mean');
    kni_corr(idx_remove,k)=nan;
    IDX_REMOVE{3,k}=idx_remove;
end
kni_corr_x=kni_corr.*L;
% plot(age,kni_corr,'.'); %inspect

kr=pp_kr';
Nb=size(kr,2);
kr_corr=nan(size(kr));
for k=1:Nb
    idx_keep=~isnan(kr(:,k));
    fit_param=polyfit(age(idx_keep),kr(idx_keep,k),1);
    kr_corr(:,k)=kr(:,k)-(fit_param(1)*(age-t0));
    [~,idx_remove]=rmoutliers(kr_corr(:,k),'mean');
    kr_corr(idx_remove,k)=nan;
    IDX_REMOVE{4,k}=idx_remove;
end
kr_corr_x=kr_corr.*L;
%% plot sigma_f and "anchoring" limit (Fig.2F);
col_hb=[0 114 189]./256;
col_gt=[217 83 25]./256;
col_kni=[0 128 0]./256;
col_kr=[126 47 142]./256;

%"Anchoring" limit
x_fine=linspace(0,480,200);
sigma_A_fine=x_fine.*std(1./L);
sigma_P_fine=sqrt((mean(L)-x_fine).^2).*std(1./L);
sigma_AP_fine=sqrt(1./(1./sigma_A_fine.^2+1./sigma_P_fine.^2));
xs_fine=x_fine./mean(L);

stats=bootstrp(1000,@(x) x_fine.*std(1./x),L);
sigma_A_fine2=nanmean(stats);
sigma_A_fine_err=nanstd(stats);

stats=bootstrp(1000,@(x) sqrt((mean(L)-x_fine).^2).*std(1./x),L);
sigma_P_fine2=nanmean(stats);
sigma_P_fine_err=nanstd(stats);

stats=bootstrp(1000,@(x) sqrt(1./(1./(x_fine.*std(1./x)).^2+1./(sqrt((mean(L)-x_fine).^2).*std(1./x)).^2)),L);
sigma_AP_fine2=nanmean(stats);
sigma_AP_fine_err=nanstd(stats);

hf=figure(10);
hold on

pgon1 = polyshape([xs_fine fliplr(xs_fine)],[sigma_A_fine+sigma_A_fine_err, fliplr(sigma_A_fine-sigma_A_fine_err)]);
hp_pgon1=plot(pgon1,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.5,'edgecolor','none');
hp4=plot(xs_fine,sigma_A_fine,'k-');

pgon2 = polyshape([xs_fine fliplr(xs_fine)],[sigma_P_fine+sigma_P_fine_err, fliplr(sigma_P_fine-sigma_P_fine_err)]);
hp_pgon2=plot(pgon2,'FaceColor',[1 0.5 0.5],'FaceAlpha',0.5,'edgecolor','none');
hp5=plot(xs_fine,sigma_P_fine,'k-');

pgon3 = polyshape([xs_fine fliplr(xs_fine)],[sigma_AP_fine+sigma_AP_fine_err, fliplr(sigma_AP_fine-sigma_AP_fine_err)]);
hp_pgon3=plot(pgon3,'FaceColor','r','FaceAlpha',0.25,'edgecolor','none');
hp6=plot(xs_fine,sigma_AP_fine,'k-');

Nbootsam=100;
marksize=3;
P_corr=hb_corr(:,1:3);
x_rel=nanmean(P_corr,1);
dx_rel=nanstd(P_corr,[],1);
stats = bootstrp(Nbootsam,@(x)nanstd(x),P_corr);
dx_rel_err=std(stats);
hp4=errorbar(x_rel,dx_rel,dx_rel_err,'o','markerfacecolor',col_hb,'color',col_hb,'markersize',marksize,'linewidth',1);

P_corr=gt_corr;
x_rel=nanmean(P_corr,1);
dx_rel=nanstd(P_corr,[],1);
stats = bootstrp(Nbootsam,@(x)nanstd(x),P_corr);
dx_rel_err=std(stats);
hp5=errorbar(x_rel,dx_rel,dx_rel_err,'s','markerfacecolor',col_gt,'color',col_gt,'markersize',marksize,'linewidth',1);

P_corr=kni_corr;
x_rel=nanmean(P_corr,1);
dx_rel=nanstd(P_corr,[],1);
[stats,bootsam] = bootstrp(Nbootsam,@(x)nanstd(x),P_corr);
dx_rel_err=std(stats);
hp6=errorbar(x_rel,dx_rel,dx_rel_err,'d','markerfacecolor',col_kni,'color',col_kni,'markersize',marksize,'linewidth',1);

P_corr=kr_corr;
x_rel=nanmean(P_corr,1);
dx_rel=nanstd(P_corr,[],1);
stats = bootstrp(Nbootsam,@(x)nanstd(x),P_corr);
dx_rel_err=std(stats);
hp7=errorbar(x_rel,dx_rel,dx_rel_err,'h','markerfacecolor',col_kr,'color',col_kr,'markersize',marksize,'linewidth',1);

plot([0 1],0.016*[1 1],'k:');
plot([0 1],0.008*[1 1],'k--');
grid on
hold off
axis([ 0 1 0 0.02]);
hl=legend([hp4, hp5, hp6, hp7],{'{hb}','{gt}','{kni}','{kr}'},'interpreter','tex','location','northeast');

xlabel('x/L');
ylabel('\sigma_{x/L}');

set(gca,'fontsize',14,'fontname','helvetica','tickdir','in','linewidth',1);
box off
pbaspect([2.5 1 1]);

