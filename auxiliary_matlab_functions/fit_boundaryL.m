function [peak_position,param] = fit_boundaryL(x, I)
% fit boundary left (low to high)
% author: mnikolic@princeton.edu

%% Fit: boundary midpoint going from low to high as a function of x
[xdata, ydata] = prepareCurveData( x, I );
%%
max_X=max(xdata);
max_Y=max(ydata);
min_Y=min(ydata);
min_x=min(xdata);

if isrow(ydata)
    ydata=ydata';
end
%%
%find the peak of the provided y data.
smooth_ydata=smooth(ydata,100);

mid_point_y=0.5*(max_Y+min_Y);
[mY,loc,w,prominence]=findpeaks(smooth_ydata,'MinPeakHeight',0.5*(max_Y+min_Y));
[~,idx]=max(prominence);
loc=loc(idx);

%%

xdata(isnan(ydata))=[];
ydata(isnan(ydata))=[];

%use data until the peak position
ydata2=ydata(1:loc);
xdata2=xdata(1:loc);

ydata2=smooth(ydata2,30);
xfine=linspace(xdata2(1),xdata2(end),length(xdata2)*100);
yfine=interp1(xdata2,ydata2,xfine,'linear');
tol=(5e-4)*(max_Y-min_Y);
loc_mid=find(abs(yfine-mid_point_y)<tol);
loc_mid=min(loc_mid);

peak_position=xfine(loc_mid);


%residual code to maintain behavior of functions that use this function. 
%Previously these were the parameters of the sigmoid fit (erf). They give
%the same results. 
param=[peak_position; xdata(loc); mid_point_y; min_Y];
%param here are not the same as the usual params of fitting the erf, but
%they are recorded in a similar way nonetheless. param(1)=peak position
% param(2) = proxy for peak width (distance from midpoint to maximum).
% param(3) is the value of gene expression (intensity) at the midpoint
% param(4) is the background but in this case it is just the minimum of the
% intensity in the profile.
end



