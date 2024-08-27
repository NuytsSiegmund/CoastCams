function [B]=ImagePreProcessing_20090121Taller(A,icmin,icmax,dt,resc,methodPreT)

% size(A) = n x m
[nt, nc, ncol] = size(A);

% Spatial resolution
resc = round(nc./100);

clear B;
A3 = A;

if(methodPreT == 0)
    B=double(A);
elseif(methodPreT==1)
    for ic=icmin:resc:icmax
        Tcoupure = 1.5; %Cut-off frequency
        fr = 1/dt;
        Val = (1/Tcoupure)*2*(1/fr);
        ord = 1000; 
        fil = fir1(ord,Val,'low');
        kk1 = conv(fil,A3(:,ic));
        S = detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));

        Tcoupure = 20; %Cut-off frequency
        fr = 1/dt;
        Val = (1/Tcoupure)*2*(1/fr);
        ord = 1000; 
        fil = fir1(ord,Val,'high');
        kk1 = conv(fil,S);
        y = detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));
     
        B(:,ic) = y;
    end

    elseif(methodPreT == 2)
        for ic=icmin:4*resc:icmax   
         [hs,htiers,hrms,trms,Hmax,h] = Wave_Char(smooth(double(A(:,ic,1)),2),dt,0,1);   
            T(ic) = trms;
        end
        To = min(T(find(T>3)))
        
              
        for ic = icmin:resc:icmax        
            v = double(A(:,ic,1));
            B(:,ic) = smooth (v-smooth(v,round(6*To)),round((2*To)/3));
        end
        
elseif(methodPreT == 3)
    for ic = icmin:resc:icmax
        v = detrend(double(A(:,ic,1)));
        B(:,ic) = smooth(smooth(v-smooth(v,60),10),3);
    end
    
elseif(methodPreT == 4)
    clear A2
    A2(1:nt,icmin:icmax)=zeroavg(double(A(1:nt,icmin:icmax)));
    for ic = icmin:icmax
        count1 = timeseries(A2(:,ic),(1:nt)*dt); %Put extrema for unfiltered signal
        interval = [0.12 0.14];
        y1 = idealfilter(count1,interval,'pass');
        y = double(y1.Data);
        B(:,ic) = y;
        B(:,ic)=Norma(B(:,ic),2);
    end
elseif(methodPreT == 5)
    for ic = icmin:resc:icmax 
        v = double(A(:,ic,1));
        B(:,ic) = smooth(v-smooth(v',7),2);
        B(:,ic) = Norma(B(:,ic),2);
    end
else
facth = 1;
factv = 1;

B = double(A);
B(factv+1:nt-factv-2,facth+1:nc-facth-2 )= smoothc(double(A),factv,facth);
end

ii=find(abs(mean(B))>0&isnan(abs(mean(B)))==0);

for irt = 1:nt
    try
     B(irt,icmin:icmax)=interp1(ii,B(irt,ii),icmin:icmax);
    end
end

disp('Pre-Processing OK')
end

function [lmval,indd]=lmax(xx,filt)
%LMAX 	[lmval, indd]=lmax(xx,filt). Find local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the 
%	corresponding indexes, FILT is the number of passes of the small
%	running average filter in order to get rid of small peaks.  Default
%	value FILT =0 (no filtering). FILT in the range from 1 to 3 is 
%	usially sufficient to remove most of a small peaks
%	For example:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	[b,a]=lmax(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMIN, MAX, MIN
	
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filt=0; 
	else
x1=x(1); x2=x(len_x); 
	for jj=1:filt,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end
lmval=[]; indd=[];
i=2;		% start at second data point in time series
    while i < len_x-1,
	if x(i) > x(i-1)
	   if x(i) > x(i+1)	% definite max
lmval =[lmval x(i)];
indd = [ indd i];
	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];  	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 2;  		% skip 2 points
	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end
if filt>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end
	  for ii=1:length(indd), 	% Find the real maximum value
	    [val(ii) iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
	    iind(ii)=indd(ii) + iind(ii)  -rng-1;
	  end
  indd=iind; lmval=val;
else
end
end

function [N,Ma]=zeroavg(M,c)

% function [N,Ma]=zeroavg(M,c)
%
% This function computes the average of each column of a matrix M
% and subtracts it % from every entry in that column. 
% If the column contains a time % serie, the time serie will become 
% zero averaged. 
%
% If an extra argument "c"
% is passed to the function 
% (independent of its value) the averageing will
% not be performed over its columns bit over its rows.
% The zero averaged matrix is returns in matrix N.

[p,q]=size(M);

if nargin == 1
   Ma=sum(M)/p;
   N = M - repmat(Ma,p,1);
end
if nargin>1
   Ma=(sum(M')/q)';
   N=M-repmat(Ma,1,q);
end
end

% Normalise un vecteur presentant 
%une variation d'amplitude par une moyenne glissante
%USAGE [f2]=Norm(f,n)
function Normf=Norma(f,n)

F=f;
version=3;
%n Valeur pour la moyenne glissante
%n=3;

%valeur absolue de f

Normf=F.*0;

if(version==1)
%
clear fct2

for it=1:length(f)
        for i=1:size(F,2)
        f=F(:,i);

    if(it<=n)
  Normft(it)=f(it)./nanmean(abs(f(it:it+n)));
elseif(it>=length(f)-n)
Normft(it)=f(it)./nanmean(abs(f(it-n:it)));
else
  Normft(it)=f(it)./nanmean(abs(f(it-n:it+n)));  
    end
Normf(1:length(Normft),i)=Normft;
        end
end
elseif(version==2)

    for i=1:size(F,2)
        f=F(:,i);
    % Version matricielle
f2=[f(1) ;f(1); f; f(length(f)); f(length(f))];   
absf=[abs(f(1)) ;abs(f(1)); abs(f); abs(f(length(f))); abs(f(length(f)))];
absp2=[abs(f(1)); abs(f(1)); abs(f(1)) ;abs(f(1)); abs(f)];
absp1=[abs(f(1)) ;abs(f(1)); abs(f(1)); abs(f); abs(f(length(f)))];
absm1=[abs(f(1)) ;abs(f) ;abs(f(length(f))) ;abs(f(length(f))); abs(f(length(f)))];
absm2=[abs(f) ;abs(f(length(f))); abs(f(length(f))); abs(f(length(f))); abs(f(length(f)))];

f2=5*f2./((absf+absp2+absp1+absm1+absm2));
Normft=f2(3:(length(f2)-2));
Normf(:,i)=Normft;
    end
    
    
elseif(version==3)
    
for ir=1:size(F,2)
t22=hilbert(F(:,ir));
Normf(:,ir)=F(:,ir)./smooth(sqrt(t22.*conj(t22)),18*2);
end
    
end %if version
end
