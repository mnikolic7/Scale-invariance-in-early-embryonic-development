%This script loads the gap gene data and removes (manually) the
%experimental sessions with high background, returning the N=301 gap gene
%profiles normalized and used for calcluation of dI. See SFig.3 in SI
%Appendix.

%author: mnikolic@princeton.edu
clear;
clc;
%% load gap gene data
load('rawProfiles_gapGenes_Hb_Gt_Kni_Kr.mat');

%%
exp_session=[data.ExpSession_idx];
idx2=find(exp_session==2);
idx4=find(exp_session==4);
idx5=find(exp_session==5);
idx8=find(exp_session==8);
data([idx2,idx4,idx5,idx8])=[];
%remove experimental sessions that have higher background (See SFig.3)
exp_session=[data.ExpSession_idx];

%
Nemb=length(data);
em_idx=[data.index];
age=[data.age];
L=[data.L];

%raw profiles.
profiles_hb_raw=horzcat(data.Hb);
profiles_gt_raw=horzcat(data.Gt);
profiles_kni_raw=horzcat(data.Kni);
profiles_kr_raw=horzcat(data.Kr);

xs=data(1).xs;
thickness=[data.dist]; %cell membrane thickness
Nexp=length(unique(exp_session));
end_idx_of_each_exp=find(diff(exp_session)>0);
end_idx_of_each_exp=[end_idx_of_each_exp,Nemb];
ES=unique(exp_session);

%%
% collected profiles after time correcting as in Dubuis+al_MSB_2013. Time
% correction is applied to all embryos within a given experimental session,
% and the (expression) profiles are then pooled in a larger set which
% amounts to N=301 total embryos. 
collected_hb=nan(900,301); 
collected_gt=nan(900,301);
collected_kni=nan(900,301);
collected_kr=nan(900,301);

% initialize new variables for corresponding pooled L, age, etc. and keep
% track of them.
collected_L=nan(301,1);
collected_age=nan(301,1);
collected_em_idx=nan(301,1);
collected_thickness=nan(301,1);

for ii=1:Nexp
    idx=find(exp_session==ES(ii));
    
    curr_L=L(idx);
    curr_age=age(idx);
    curr_thickness=thickness(idx);
    curr_em_idx=em_idx(idx);
    
    curr_hb=profiles_hb_raw(:,idx);
    curr_gt=profiles_gt_raw(:,idx);
    curr_kni=profiles_kni_raw(:,idx);
    curr_kr=profiles_kr_raw(:,idx);
    
    [profiles_kni3,curr_thickness,curr_L,curr_age,curr_em_idx]=correct_time(curr_kni,curr_thickness,curr_L,curr_age,curr_em_idx);
    [profiles_kr3,curr_thickness,curr_L,curr_age,curr_em_idx]=correct_time(curr_kr,curr_thickness,curr_L,curr_age,curr_em_idx);
    [profiles_gt3,curr_thickness,curr_L,curr_age,curr_em_idx]=correct_time(curr_gt,curr_thickness,curr_L,curr_age,curr_em_idx);
    [profiles_hb3,curr_thickness,curr_L,curr_age,curr_em_idx]=correct_time(curr_hb,curr_thickness,curr_L,curr_age,curr_em_idx);
    
    [profiles_kni3,mean_kni,std_kni]=normalize_profiles(profiles_kni3);
    [profiles_kr3,mean_kr,std_kr]=normalize_profiles(profiles_kr3);
    [profiles_hb3,mean_hb,std_hb]=normalize_profiles(profiles_hb3);
    [profiles_gt3,mean_gt,std_gt]=normalize_profiles(profiles_gt3);
    
    % record the mean profile and the stdev of gene expression as a
    % function of xs in each experimental session (ii)
    mean_profiles_hb(:,ii)=nanmean(profiles_hb3,2);
    std_profiles_hb(:,ii)=nanstd(profiles_hb3,[],2);
    
    mean_profiles_gt(:,ii)=nanmean(profiles_gt3,2);
    std_profiles_gt(:,ii)=nanstd(profiles_gt3,[],2);
    
    mean_profiles_kni(:,ii)=nanmean(profiles_kni3,2);
    std_profiles_kni(:,ii)=nanstd(profiles_kni3,[],2);
    
    mean_profiles_kr(:,ii)=nanmean(profiles_kr3,2);
    std_profiles_kr(:,ii)=nanstd(profiles_kr3,[],2);
    
    collected_hb(:,idx)=profiles_hb3;
    collected_gt(:,idx)=profiles_gt3;
    collected_kni(:,idx)=profiles_kni3;
    collected_kr(:,idx)=profiles_kr3;
    
    collected_L(idx)=curr_L;
    collected_age(idx)=curr_age;
    collected_em_idx(idx)=curr_em_idx;
    collected_thickness(idx)=curr_thickness;
end
%% plot variance of intensity divided by total amplitude squared. SFig.3
cmap=turbo;
def_colors=cmap(round(128:(128/4):256),:);

hf=figure(1);

mean_profiles=mean_profiles_hb;
var_profiles=std_profiles_hb.^2;
maxI=max(mean_profiles,[],1);
minI=min(mean_profiles,[],1);
subplot(2,2,1)
hold on
for k=1:Nexp
        c='b';
        plot(xs,var_profiles(:,k)./(maxI(k)-minI(k)).^2,'-','color',c)
end
plot([0 1],0.001*[1 1],'k--');
plot([0 1],0.1*[1 1],'k-.');
hold off
xlabel('x/L')
ylabel('{\boldmath$\sigma^2_g=\sigma^2_I/(I_{max}- I_{min})^2$}','interpreter','latex')
set(gca,'yscale','log')
ylim(10.^[-5 -1])
set(gca,'fontname','helvetica','fontsize',12,'linewidth',1.25)
title('Hb')

mean_profiles=mean_profiles_hb;
var_profiles=std_profiles_hb.^2;
maxI=max(mean_profiles,[],1);
minI=min(mean_profiles,[],1);
subplot(2,2,1)
hold on
for k=1:Nexp
        c='b';
        plot(xs,var_profiles(:,k)./(maxI(k)-minI(k)).^2,'-','color',c)
end
plot([0 1],0.001*[1 1],'k--');
plot([0 1],0.1*[1 1],'k-.');
hold off
xlabel('x/L')
ylabel('{\boldmath$\sigma^2_g=\sigma^2_I/(I_{max}- I_{min})^2$}','interpreter','latex')
set(gca,'yscale','log')
ylim(10.^[-5 -1])
set(gca,'fontname','helvetica','fontsize',12,'linewidth',1.25)
title('Hb')

mean_profiles=mean_profiles_gt;
var_profiles=std_profiles_gt.^2;
maxI=max(mean_profiles,[],1);
minI=min(mean_profiles,[],1);
subplot(2,2,2)
hold on
for k=1:Nexp
        c='b';
        plot(xs,var_profiles(:,k)./(maxI(k)-minI(k)).^2,'-','color',c)
end
plot([0 1],0.001*[1 1],'k--');
plot([0 1],0.1*[1 1],'k-.');
hold off
xlabel('x/L')
ylabel('{\boldmath$\sigma^2_g=\sigma^2_I/(I_{max}- I_{min})^2$}','interpreter','latex')
set(gca,'yscale','log')
ylim(10.^[-5 -1])
set(gca,'fontname','helvetica','fontsize',12,'linewidth',1.25)
title('Gt')

mean_profiles=mean_profiles_kni;
var_profiles=std_profiles_kni.^2;
maxI=max(mean_profiles,[],1);
minI=min(mean_profiles,[],1);
subplot(2,2,3)
hold on
for k=1:Nexp
        c='b';
        plot(xs,var_profiles(:,k)./(maxI(k)-minI(k)).^2,'-','color',c)
end
plot([0 1],0.001*[1 1],'k--');
plot([0 1],0.1*[1 1],'k-.');
hold off
xlabel('x/L')
ylabel('{\boldmath$\sigma^2_g=\sigma^2_I/(I_{max}- I_{min})^2$}','interpreter','latex')
set(gca,'yscale','log')
ylim(10.^[-5 -1])
set(gca,'fontname','helvetica','fontsize',12,'linewidth',1.25)
title('Kni')

mean_profiles=mean_profiles_kr;
var_profiles=std_profiles_kr.^2;
maxI=max(mean_profiles,[],1);
minI=min(mean_profiles,[],1);
subplot(2,2,4)
hold on
for k=1:Nexp
        c='b';
        plot(xs,var_profiles(:,k)./(maxI(k)-minI(k)).^2,'-','color',c)
end
plot([0 1],0.001*[1 1],'k--');
plot([0 1],0.1*[1 1],'k-.');
hold off
xlabel('x/L')
ylabel('{\boldmath$\sigma^2_g=\sigma^2_I/(I_{max}- I_{min})^2$}','interpreter','latex')
set(gca,'yscale','log')
ylim(10.^[-5 -1])
set(gca,'fontname','helvetica','fontsize',14,'linewidth',1.25)
title('Kr')

%% prepare data for calculating information gain.


%middle part of Kni from 0.15 to 0.45 includes the small middle stripe of
%Kni which varies systematically along the dorsoventral axis. Remove this
%part from the analysis, as we cannot exclude that the slight error in DV
%orientation of the embryo during imaging does not correlate with length,
%which would affect the entropy calculation. 
% See SI Appendix, section D. See also
% https://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=1&ftext=FBgn0001320
% for an example of the Kni expression profile.
collected_kni(100:400,:)=nan; %note that the indices 100-400 correspond to xs from 0.151-0.45 (because index 1 = 0.051)



%divide profiles into position bins of witdth 0.01L
profiles_kr100=nan(90,Nemb);
profiles_kni100=nan(90,Nemb);
profiles_gt100=nan(90,Nemb);
profiles_hb100=nan(90,Nemb);

for j=1:301
    for k=1:90
        idx_sel=10*k-5;
        idx_mean=(10*k)+(-9:1:0);
        profiles_kr100(k,j)=nanmean(collected_kr(idx_mean,j) ,1);
        profiles_kni100(k,j)=nanmean(collected_kni(idx_mean,j) ,1);
        profiles_gt100(k,j)=nanmean(collected_gt(idx_mean,j) ,1);
        profiles_hb100(k,j)=nanmean(collected_hb(idx_mean,j) ,1);
        xs100(k)=mean(xs(idx_mean));
    end
end

%inspect
figure(2)
subplot(2,2,1)
plot(xs100,profiles_hb100,'.-')
xlabel('xs');
ylabel('g_{\rm Hb}');
subplot(2,2,2)
plot(xs100,profiles_gt100,'.-')
xlabel('xs');
ylabel('g_{\rm Gt}');
subplot(2,2,3)
plot(xs100,profiles_kni100,'.-')
xlabel('xs');
ylabel('g_{\rm Kni}');
subplot(2,2,4)
plot(xs100,profiles_kr100,'.-')
xlabel('xs');
ylabel('g_{\rm Kr}');
profiles=profiles_kni100;
%%
data=struct('age',[],'L',[],'Hb',[],'Gt',[],'Kni',[],'Kr',[]);
for k=1:301
    data(k).age=collected_age(k);
    data(k).L=collected_L(k);
    data(k).Hb=collected_hb(:,k);
    data(k).Gt=collected_gt(:,k);
    data(k).Kni=collected_kni(:,k);
    data(k).Kr=collected_kr(:,k);
end