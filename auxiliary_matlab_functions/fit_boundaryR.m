function [peak_position,param] = fit_boundaryR(x, I)
% fit boundary right (high to low)
% author: mnikolic@princeton.edu

%% Fit: 
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

xdata(isnan(ydata))=[];
ydata(isnan(ydata))=[];

ydata2=ydata;
xdata2=xdata;
mid_point_y=0.5*(max_Y+min_Y);

ydata2=smooth(ydata2,30);

xfine=linspace(xdata2(1),xdata2(end),length(xdata2)*100);
yfine=interp1(xdata2,ydata2,xfine,'spline');
tol=(1e-3)*(max_Y-min_Y);
loc_mid=find(abs(yfine-mid_point_y)<tol);
loc_mid=max(loc_mid);

peak_position=xfine(loc_mid);

%residual code to maintain behavior of functions that use this function. 
%Previously these were the parameters of the sigmoid fit (erf). They give
%the same results. 
param=[peak_position; peak_position; mid_point_y; min_Y];
%param here are not the same as the usual params of fitting the erf, but
%they are recorded in a similar way nonetheless. param(1)=peak position
% param(2) = proxy for peak width (distance from midpoint to maximum).
% param(3) is the value of gene expression (intensity) at the midpoint
% param(4) is the background but in this case it is just the minimum of the
% intensity in the profile.
end



