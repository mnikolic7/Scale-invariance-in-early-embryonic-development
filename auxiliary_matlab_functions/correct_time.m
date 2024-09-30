function [profiles_corrected,t_cell,L,age,em_idx] = correct_time(profiles,t_cell,L,age,em_idx)
%following Dubuis et al. 2013, Mol. Sys. Bio. 
% This function correct the timing of
% gap gene profiles measured in a given experiment.
% profiles is a matrix with 1000 rows (one for each x/L AP position) and
% each column is a separate embryo. There are Nemb embryos (columns). 
% t_cell is a vector of length Nemb that contains cellularization distances
% in microns. Output is a matrix of corrected profiles, same size and shape
% as the provided variable profiles.
% returned profiles are sorted by t_cell. For this reason, you should
% provide L, and age so that they can be sorted, all in the same place.  
%   author: mnikolic@princeton.edu

% %%
% profiles=profiles_gt2;
% t_cell=thickness_clean;
% L=L_clean;
% age=age_clean;
%% sort by thickness.
Nemb=size(profiles,2);
Nx=size(profiles,1);
profiles_sorted=nan(Nx,Nemb);
profiles_corrected=nan(Nx,Nemb);
profiles_trend=nan(Nx,Nemb);
[t_cell,idx_sort]=sort(t_cell,'ascend');
%keep track of the corresponding auxiliary variables L, age, and embryo index.
L=L(idx_sort);
em_idx=em_idx(idx_sort);
age=age(idx_sort);
for ii=1:Nx
    profiles_sorted(ii,:)=profiles(ii,idx_sort);
end

%% correct for time
filt_wid=2*Nemb; %Width of the (weighted!) filter. This number does not matter much. The next line determines the effective filt width. 2*Nemb ensures all embryos are included.
filt_sigma=2.5; %in microns for thickness. (t_cell is in microns). %as in Dubuis+al_2013.
for ii=1:Nx %for each position x

    for jj=1:Nemb%/4*3 %go through each position for smoothing. 
        idx_filt=jj-round(filt_wid./2):jj+round(filt_wid/2); % allows for custom filter width, but in here it is set to be 2x the length of the embryo - meaning for each position it is always contains the entire array.
        idx_filt(idx_filt<=0)=[]; %remove impossible indices <=0 and >Nemb
        idx_filt(idx_filt>Nemb)=[];
        

        center=t_cell(jj);
        t_filt=t_cell(idx_filt);
        I_filt=profiles_sorted(ii,idx_filt);
        weights=exp(-(t_filt-center).^2./(2*filt_sigma.^2));
        weights=weights./sum(weights);
        if isrow(weights)
            weights=weights';
        end
        corr_intensity=I_filt*weights; %this is dot product.   
        profiles_trend(ii,jj)=corr_intensity;

    end
end

%% correct only for the temporal trend
t0=45;
idx_t0=find(min(abs(age-t0)));
profiles_at_t0=repmat(profiles_trend(:,idx_t0),1,Nemb);

%subtract only the deviation of the trend from the t0 reference, leaving
%the biological variation in. See Dubuis+al_2013 for more details.
profiles_corrected=profiles_sorted-(profiles_trend-profiles_at_t0);

end