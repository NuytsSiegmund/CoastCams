function  fsig=FilterMean(a,b)
% Filter data using running mean filter
% a: input vector
% b: filter width

version = 1;
[n1, n2] = size(a);

if n2 == 1
    a = a';
end

% Version 1
if (version == 1)
    la = length(a);

    if (la <(2 * b))
        disp('Attention: vector too small for filtering')
    end
    
    fsig(1:la) = 0;
    fsig2(1:la + 2 * b) = 0;
    
    fsig2(1:b) = mean(a(1:b));
    fsig2(b + 1:b + la) = a;
    fsig2(b + la + 1:length(fsig2)) = mean(a(la - b:la));
    
    k = 0;
    
    for kk = b + 1:b + la
        k = k + 1;
        fsig(k) = mean(fsig2(kk - b:kk + b));
    end
    
    % Uncomment the following two lines to add additional filtering at the edges
    % fsig(la - (b - 1):la) = (a(la - 1) + a(la - (b - 1))) / 2;
    % fsig(1:(b - 1)) = (a(1) + a((b - 1))) / 2;
elseif(version==2)
    
    a=double(a);
    
    a2=a';
    
am(1:b)=a(1);
ap(1:b)=a(length(a));
a2=[am' ;a'; ap'];
a3(1:length(a))=a';
c=0;
for i=1:b
a3=a3+a2((b+1:(length(a2))-b)+i)'+a2((b+1:(length(a2))-b)-i)';
c=c+1;
end
fsig=a3/(2*b+1);

end % Version

fsig(find(isnan(fsig)==1))=a(find(isnan(fsig)==1));
end   