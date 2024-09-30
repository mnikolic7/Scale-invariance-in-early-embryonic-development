function [peak_position,param] = fit_peak_splitLRmidpoint(x, I)
% identify the peak, then split it in two halfs and find the two intensity
% midpoints. Return two positions: left and right midpoint of the larger
% (peak shaped) domain. Used for Kruppel expression profile, for example. 
% author: mnikolic@princeton.edu

%%
[xdata, ydata] = prepareCurveData( x, I );

%%
max_X=max(xdata);
max_Y=max(ydata);
min_Y=min(ydata);
min_X=min(xdata);

if isrow(ydata)
    ydata=ydata';
end
if isrow(xdata)
    xdata=xdata';
end
%% find large peak
smooth_ydata=smooth(ydata,100);

mid_point_y=0.5*(max_Y+min_Y);
[mY,loc,w,prominence]=findpeaks(smooth_ydata,'MinPeakHeight',0.5*(max_Y+min_Y));
[~,idx]=max(prominence);
loc=loc(idx);
%inspect
% plot(xdata,smooth_ydata,'.-');
% hold on
% plot(xdata,ydata,'.-')
% plot(xdata,ones(size(xdata))*min_Y,'k-');
% plot(xdata,ones(size(xdata))*max_Y,'k-');
% plot(xdata,ones(size(xdata))*(max_Y+min_Y)*0.5,'k-');
% plot(xdata(loc)*[ 1 1], [0 1.3],'k-');
% hold off

center=xdata(loc);
xdataL=xdata(1:loc);
ydataL=ydata(1:loc);

xdataR=xdata(loc+1:end);
ydataR=ydata(loc+1:end);

mid_point_y=0.5*(max_Y+min_Y);
%% find midpoints 
xfine=linspace(xdataL(1),xdataL(end),length(xdataL)*100);
yfine=interp1(xdataL,ydataL,xfine,'spline');
tol=(1e-3)*(max_Y-min_Y);
loc_mid=find(abs(yfine-mid_point_y)<tol); %find where the intensity=1/2 maximum.
if isempty(loc_mid) %if you can't find mid point, extrapolate a line on the 30% of the interval and find where it would hit the middle. 
    npts=length(xdataL)
    idx_pts=round((0.65*npts):(0.85*npts));
    options = fitoptions('Method','NonlinearLeastSquares', 'Robust', 'LAR','StartPoint',[min_Y,2*(max_Y-min_Y)/(max_X-min_X)]);
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xdataL(idx_pts), ydataL(idx_pts), ft, options);
    x_mp=(mid_point_y-fitresult.p2)./fitresult.p1
    if x_mp<=min_X
        x_mp=nan;
    end
else
    %
    loc_mid=round(max(loc_mid));
    x_mp=xfine(loc_mid);
end
peak_position=x_mp;

paramL=[peak_position; abs(center-peak_position); mid_point_y; min_Y];
%%
xfine=linspace(xdataR(1),xdataR(end),length(xdataR)*100);
yfine=interp1(xdataR,ydataR,xfine,'spline');
tol=(1e-3)*(max_Y-min_Y);
loc_mid=find(abs(yfine-mid_point_y)<tol);
loc_mid=round(min(loc_mid));

peak_position=xfine(loc_mid);

paramR=[peak_position; abs(center-peak_position); mid_point_y; min_Y];

%%
param_final=[paramL;paramR];
peak_position=param_final(1:4:end);
param=param_final;
end