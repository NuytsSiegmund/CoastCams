function [hs,htiers,hrms,trms,Hmax,h,Tp,tindiv]=Wave_Char(d,dt,filt,meth)
% Calculate wave characteristics from a time series of water level
% Inputs
% d: time series vector
% dt: time step
% filt: filter option (1 for filtering, 0 for no filtering)
% meth: Mean zero crossing option
%       meth=1 mean zero crossing
%       meth=2 up zero crossing
%       meth=3 down zero crossing

% Outputs
% hs: significant wave height
% htiers: wave height of the upper third of the waves
% hrms: root mean square wave height
% trms: root mean square wave period
% Hmax: maximum wave height
% h: individual wave heights (for example, to make a histogram)
% Tp: peak wave period
% tindiv: individual wave periods (for example, to make a histogram)
 
% Number of measurement time steps (nt)
[n1, n2] = size(d);
if n1 == 1
    d = d';
end
sz = size(d);
nt = sz(1);
zc = meth;

% Detrend the raw signal
k = d(1:nt);

x = 1:nt;
a = reshape(x, nt, 1);
p = polyfit(a, k, 1);
k2 = k - p(1).*a;
km = k2;
k2 = km - mean(km);

clear k4 
% Apply a filter if needed
if (filt == 1)
    freq = (1./dt)*(1:nt)/nt;
    fcb = 0.01;
    fch = 1/0.1;
    f1 = fft(k2);
    Indb = find(freq < fcb);
    Indh = find(freq > fch);
    f1(Indb) = 0.0;
    Sinv = ifft(f1);
    k3 = Sinv.*conj(Sinv);
    f1(Indh) = 0.0;
    Sinv = ifft(f1);
    k4 = real(Sinv);
    k2 = k4;
    k3 = FilterMean(k2, round(dt));
else
    k3 = FilterMean(k2, round(dt));
    k4 = k2;
end

% Localize the zeros of the signal
j = 0;
z = zeros(1, floor(nt./2));  % Pre-allocate z with zeros

if (zc == 1)
    % Mean zero crossing
    for i = 2:nt-1
        if (k3(i)*k3(i+1) <= 0 && k3(i-1)*k3(i) >= 0)
            j = j+1;
            z(j) = i;
        end
    end
elseif (zc == 2)
    % Up zero crossing
    for i = 2:nt-1
        if (k3(i)*k3(i+1) <= 0 && mean([k3(i-1) k3(i)]) < 0)
            j = j+1;
            z(j) = i;
        end
    end   
elseif (zc == 3)
    % Down zero crossing
    for i=2:nt-1
        if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])>0 ) 
            j=j+1;
            z(j)=i;
        end
    end    
end

% Remove unused elements from z
z = z(1:j);

% Calculate the size of the vector that stores the indices of zeros
k = 0;
sz = size(z);
for i = 1:sz(2)
    if (z(i)~=0)
    k = k + 1;
    end
end

% Find the max height of the wave between two zeros
hmax(1:k-1) = 0;
hmin(1:k-1) = 0;

for i=1:k-1
    if(mean(k4(z(i):z(i+1)))>0.)
        if(max(k3(z(i):z(i+1)))>0)
    hmax(i)=max(k4(z(i):z(i+1)));
        else
    hmax(i)=0;
        end
    hmin(i)=0.;
    
else
    hmax(i)=0.;
        if(min(k3(z(i):z(i+1)))<0)
        hmin(i)=min(k4(z(i):z(i+1)));
    
        else
        hmin(i)=0.; 
        end    
    end
end

% Calculate the height of each wave passing through the sensor
h = [];
kh = 0;
h(1:k-2) = 0;

for i=1:k-2
    
 if(hmin(i)~=0 & hmax(i+1)~=0)
        if(hmin(i+1)==0 & hmax(i)==0)        
        kh=kh+1;
     h(kh)=hmax(i+1)-hmin(i);
        end
 end
end

% Calculate the average rms height of the waves over the measurement period
clear hrms ord htiers

[sd ord]=sort(h(find(h~=0)),'descend');

%Max Waveheight
Hmax = max(h);

%1/3 Waveheight
htiers=mean(h(ord(1:round(length(ord)/3))));

%RMS Waveheight
hrms=mean(h(1:kh));

%Significant Waveheight
hs=4.*std(k4);

clear t
t(1:k-1)=0;
for i=1:k-1
    if(z(i)~=0 | z(i+1)~=0 ) 
        t(i)=2*(z(i+1)-z(i))*dt;
    end
end

tindiv = t;

%Calculate the rms period of the waves and display it
if(zc==1)
trms=mean(t);
elseif(zc==2)
trms=mean(t)/2;
elseif(zc==3)
trms=mean(t)/2;
end    

% Peak wave period
y=detrend(d);

nx=length(y);
%Temporal resolution (sampling period) (in s)

%FFT specifying the vector size
z=fft(y);
[n1 n2]=size(z);
if n1==1
    z=z';
end

t=1./((1:nx/2)./(dt*nx));
in=find(t<6|t>14);
x=1:nx/2;
k=z(x).*conj(z(x)/nx);
k(in)=0;

[sa ind]=lmax(FilterMean(FilterMean(k,3),5));

Tp=1./(ind./(dt*nx));
end