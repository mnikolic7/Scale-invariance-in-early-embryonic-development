function [dI_at_inf,err,NN,dI,AAdI2] = calculate_dI(profiles,L,varargin)
%Input: profiles: intensity of gene expression in relative coordinates.
%       L: embryo length (corresponding to columns of profiles)
%It calculates the positional information gain from knowing the embryo
%length according to the Eq. [14] in the main text and Eq.[S22] in the 
%SI Appendix, while following the procedure in SI Appendix, section E:
%Estimating \Delta I from limited data.
%
%Output: dI_at_inf: unbiased estimated of \Delta I
%        err: error bar of dI
%        NN: number of sampled embryos in each "experiment/draw"
%        dI: calculated dI as a function of NN, each experiment/draw is 
%            repeated 100 times for each NN(i). 
%        AAdI2: quadratic fit of the systematic behavior of mean(dI)_across
%        experiments/draws.
%authors: William Bialek, 
%         Milos Nikolic (mnikolic@princeton.edu)

%%
total_N_embryos=size(profiles,2);
Nxpoints=size(profiles,1);
profiles=profiles'; %from now on, the columns are positions (row vectors) and rows are one for each embryo.
Nem = total_N_embryos;
%% If this option is on (by providing an additonal argument which is a string 'shuffleL'
% the indices of L will be permuted, so that L becomes a variable that by
% definition provides zero additional information. Used for control tests.
permuteL=0;
if nargin>2
    if strcmpi(varargin{1},'shuffleL')
        permuteL=1;
    end
end
%% check positions are in 90 bins (of width 0.01L), and put L into nL bins
nL=5; %use 5 length bins.

nXpts=1; %change this from 1 to 10 (this averages nXpts to reduce xs bin size to 0.01 (or arbitrary number)).
apaxis = [1:Nxpoints];
qx = ceil(apaxis/nXpts); %qx is bin label for position;
n_x_axis=length(unique(qx));

%split L in nL bins of adaptive size, with constant number of embryos in
%each bin.
sortL = sort(L,'ascend');
cutsL = sortL(ceil([1:nL]*Nem/nL));
qL = ones(Nem,1);
for n=2:5
    idx = find(L> cutsL(n-1) & L<= cutsL(n));
    qL(idx) = n;
end
%qL is label indicating to which bin of L current embryo belongs
if permuteL
    qL=qL(randperm(length(qL))); % shuffle L labels for controls
end

for n=1:nL
    Laxis(n) = mean(L(qL==n));
end
%%


tic
NN = [50:25:Nem floor(Nem/2) Nem]; %number of embryos used for calculation.
length_NN=length(NN);
parfor kk = 1:100 %repeat the experiment of drawing NN (variable number NN) embryos 100 times
    idx = randperm(Nem); 
    for n=1:length_NN %n is index to cycle through various sample sizes NN
        gg = profiles(idx(1:NN(n)),:); %gg is selection of random NN(n) profiles (gene expression values)
        LL = qL(idx(1:NN(n))); %LL is the list of L bins for each embryo profile in gg
        sigma2 = zeros(n_x_axis,nL); %rows=position bins, columns=L bins;
        sigma2x=zeros(1,n_x_axis); %make sure that the indices are right, if they are wrong, this will throw an out of bounds error.
        for i=1:n_x_axis %i=position
            for q = 1:nL %q=length bin
                sigma2(i,q) = var(reshape(gg(LL==q,qx==i),1,[])); %sigma2 = variance of all points in the bin L (qL) and bin x (qx)
            end
            sigma2x(i) = var(reshape(gg(:,qx==i),1,[])); %sigma2x is variance for position bin x, but for all lengths: all embryos.
        end
        sigma2 = sigma2./(sigma2x'*ones(1,nL)); % find the ratio of sigma|L vs sigma_allL
        dI(kk,n) = -0.5*nanmean(nanmean(log2(sigma2 + eps))); %difference of two entropies, averaged as in Eq.[14] and Eq.[S22]
        % kk is just one of the 100 repeats of this 'experiment' (in which we
        % randomly choose NN(n) embryos).
        % n is the index that indexes NN so in the end we get a table of 
        % dI (delta information) vs. sample size. 
    end
end
%
%difference in information decreases as a function of number of embryos,
% so AAdI_samp is the fitted slope of the 1./NN and dI for sampling
% experiment/draw
AAdI2_samp=nan(100,3);
for kk=1:size(dI,1)
    AAdI2_samp(kk,:)= polyfit(1./NN,dI(kk,:),2);
end


% AAdI = polyfit(1./NN,mean(dI),1); %linear fit
AAdI2= polyfit(1./NN,mean(dI),2); %quadratic fit.
toc   
%
err=std(dI); %errorbar from 100 'experiments'/draws
err=err(end-1)/sqrt(2); %errorbar from Nemb/2 experiment, scaled by sqrt(2).
extrap=AAdI2(:,3); %extrapolated value at 1/NN -> 0
dI_at_inf=extrap;


end