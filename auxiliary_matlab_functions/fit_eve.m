function [peak_position,param,fval,x_found_peaks] = fit_eve(y, I, varargin)
% fit eve profile with seven peaks
% you can provide param0 (for example from a previous fit) as an optional
% argument at the end. 
% author: mnikolic@princeton.edu

%% Fit: eve as seven gaussians
% [xData, yData] = prepareCurveData( y, I );
ydata=I;
xdata=y';

max_X=max(xdata);
max_Y=max(ydata);

%initial position guesses (from mean profile)
x00=[0.352,0.431,0.503,0.557,0.615,0.675,0.757]*max_X;
y00=[1 0.86 0.921 0.768 0.647 0.804 1]*max_Y;
w00=ones(7,1)*0.01;
param0=ones(7*3,1);

%find peaks, preliminary
I_curr=smoothdata(ydata,'movmean',20);
[peaks,loc]=findpeaks((I_curr),'MinPeakHeight',0.1*max(I_curr),'minpeakdistance',0.001*length(I_curr),'minpeakprominence',0.01*max(I_curr));
x_found_peaks=xdata(loc);

x0=x00;
y0=y00;
%for each found peak update the locations and intensities. Use these as the
%initial guess.
for k=1:length(loc)
    x_curr=xdata(loc(k));
    [~,idx]=min(abs(x_curr-x00));
    x0(idx)=x_curr;
    y0(idx)=peaks(k);
end


for k=1:7
    param0(3*k-2)=x0(k);
    param0(3*k-1)=w00(k);
    param0(3*k)=y0(k);
end

if nargin>=3
    param0=varargin{1};
end


%curve fit
[param,fval]=fitPeaksGauss(xdata,ydata,7,'provide',param0,'displayon');


peak_positions=param(1:3:end);
peak_position=max(peak_positions);

end

