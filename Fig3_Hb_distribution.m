% Reconstruct Fig3. from the main text
%author: Milos Nikolic, mnikolic@princeton.edu

clear;
clc;
%% load time corrected data,
load('gapGenes_timeCorrected_N301.mat');

%%
Nemb=length(data);
xs=(51:950)./1000;
L=[data.L];
age=[data.age];
profiles_hb=horzcat(data.Hb);
profiles_gt=horzcat(data.Gt);
profiles_kni=horzcat(data.Kni);
profiles_kr=horzcat(data.Kr);

%% Plot probablity density Fig3A

%%% code block needed for Fig3B.
cmapB=colormap('jet');
cmapB=cmapB(1:50:end,:);
%%%

figure(1)
profiles=profiles_hb;
dInt=0.05;
RInt=2;
Int=0:dInt:RInt;
profile_density=zeros(length(xs),length(Int));
profile_density_raw=zeros(length(xs),length(Int));

for ii=1:length(xs)
    for jj=1:Nemb
        curr_pos=xs(ii);
        curr_int=profiles(ii,jj);
        [~,idx_int_curr]=min(abs(curr_int-Int));
        profile_density(ii,idx_int_curr)=profile_density(ii,idx_int_curr)+1;
        profile_density_raw(ii,idx_int_curr)=profile_density_raw(ii,idx_int_curr)+1;
    end
    profile_density(ii,:)=profile_density(ii,:)./sum(profile_density(ii,:));
end
%
cmap=colormap('bone');
cmap=flipud(cmap);
profile_density2=profile_density./dInt;
profile_density3=profile_density2;
% profile_density3=imgaussfilt(profile_density3,0.2)                 
profile_density3(profile_density3==0)=nan;
pcolor(xs,Int,profile_density2')
% pcolor(y,Int,log10(profile_density3'));
view([0 90])
colormap(gca,cmap);
shading interp
hold on
plot(xs,mean(profiles,2),'-');
plot(0.47*[1 1],[0 2],'k--')
hold off
ylabel('g_{Hb}')
hcb=colorbar('location','westoutside');
% caxis([0 0.01])
ylabel(hcb,'Probability','fontname','helvetica','fontsize',12);
pbaspect([2,1,1])
set(gca,'fontname','helvetica','fontsize',12)
axis([0 1 0 2])
box off
% xlim([0.05 0.95])
% set(gca,'xtick',0.05:0.2:0.95)
%% P(g|x/L and L) vs P(g|x/L);

%create nL bins of length with roughly equal number of embryos in each.
nL=5; %n length bins;
sortL = sort(L,'ascend');
cutsL = sortL(ceil([1:nL]*Nemb/nL));
qL = ones(Nemb,1);
for n=2:5
    idx = find(L> cutsL(n-1) & L<= cutsL(n));
    qL(idx) = n;
end

%plot distribution at 0.475 \pm 0.01 x/L
dx=10;
pos_idx=425;
pos_indices=(pos_idx-ceil(dx/2)):((pos_idx+floor(dx/2)));
%
gg_curr=reshape(profiles(pos_indices,:),[],1);
Mgg=nanmean(gg_curr);
step_size=0.07; %in g
pd = fitdist(gg_curr,'Kernel','width',step_size); %use KDE to estimate the distribution
int_values = (0:step_size:2);
val_idx_mid=ceil(length(int_values)/2);
prob = pdf(pd,int_values); 

prob_L=nan(length(prob),nL);
for k=1:nL
    gg_curr=reshape(profiles(pos_indices,qL==k),[],1);
    Mgg=nanmean(gg_curr);
    step_size=0.07;
    pd = fitdist(gg_curr,'Kernel','width',step_size);
    int_values = (0:step_size:2);
    val_idx_mid=ceil(length(int_values)/2);
    prob_curr = pdf(pd,int_values);
    prob_L(:,k)=prob_curr;
end

dintv=step_size/2;



figure(2)
axis tight
for k=1:nL
    subplot(1,5,k)
    plot(prob,int_values,'k--')
    hold on
    pbaspect([1,2.5,1])
    plot(prob_L(:,k),int_values,'-','color',cmapB(k,:));
    title(['L=',num2str(mean(L(qL==k)),'%.0f')]);
    hold off
    xlim([0 2.2])
    ylim([0 1])
    if k>1
    set(gca,'ytick',[])
    
    end
    if k==3
        xlabel('Probability density');
    end
    if k==1
        ylabel('g_{Hb}')
    end
end
