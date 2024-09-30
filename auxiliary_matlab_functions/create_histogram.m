function [yCounts,yKDE,bin_centers] = create_histogram(data_points,bin_centers)
%CREATE_HISTOGRAM This wrapper function takes in a vector of data 
% creates a histogram and returns the values binned in provided
%bins, so that you can plot easily after. It also returns the Kernel
%density estimate with the gaussian kernel with std = width of first bin.
% both yKDE and yCounts are returned as probability (not probability density). 
%   author: mnikolic@princeton.edu
bin_width=bin_centers(2)-bin_centers(1); 
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];


% shift_data=shift_data(shift_data>6.09);
[yCounts,~]=histcounts(data_points,bin_edges,'Normalization','probability');
pdKDE = fitdist(data_points,'kernel','BandWidth',bin_width);
yKDE=pdf(pdKDE, bin_centers)*bin_width;

end

