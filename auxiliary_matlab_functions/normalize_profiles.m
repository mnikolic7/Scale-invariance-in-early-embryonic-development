function [norm_profiles,mean_profile,std_profile] = normalize_profiles(profiles)
% This function follows Petkova Cell, 2019 and normalizes the gap gene
% expression profiles in a given experiment by calculating the mean profile
% across all embryos (columns) and normalizes all profiles by subtracting
% the minimum of the mean profile and dividing by max-min of the mean
% profile. It also returns the mean_profile and the std_profiles in
% absolute (provided) units. 
% the input should be a column matrix, where each column is a gene
% expression profile. All profiles should be from the same experimental
% session.
%   author: mnikolic@princeton.edu

mean_profile=nanmean(profiles,2);
std_profile=nanstd(profiles,[],2);
minI=min(mean_profile);
maxI=max(mean_profile);

norm_profiles=(profiles-minI)./(maxI-minI);
end