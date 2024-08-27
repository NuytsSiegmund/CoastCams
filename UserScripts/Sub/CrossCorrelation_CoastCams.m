function [R2M, L2M, T2M, Hs, RM] = CrossCorrelation_Coastcams(A2, dpha, dt, dc)
% Copyright 2017 Rafael Almar (IRD, France)- rafael.almar@ird.fr
% Modified to include error checking and debugging

% Input validation
if isempty(A2) || ~isnumeric(A2)
    error('Input A2 must be a non-empty numeric matrix');
end
if ~isnumeric(dpha) || ~isscalar(dpha) || dpha <= 0
    error('dpha must be a positive scalar');
end
if ~isnumeric(dt) || ~isscalar(dt) || dt <= 0
    error('dt must be a positive scalar');
end
if ~isnumeric(dc) || ~isscalar(dc) || dc <= 0
    error('dc must be a positive scalar');
end

% Initialize outputs
R2M = []; L2M = []; T2M = []; Hs = []; RM = [];

try
    A2 = double(A2(:,:,1));
    A2 = detrend(A2);
    dc = round(dc/2) * 2;
    res = 1; % Resolution for cross-correlation
    [nt, nc, ~] = size(A2);
    n = floor(dpha/dt);

    % Check if n is valid
    if n <= 0 || n >= nt
        error('Invalid n value calculated from dpha and dt');
    end

    % Preallocate R2M
    R2M = zeros(nc - dc + 1, dc - 1);

    for ic = 1 + dc/2 : res : nc - dc/2
        R2 = zeros(1, dc - 1);
        for lc = 1 : dc - 1
            % Check for valid indices
            if ic - dc/2 > 0 && ic - dc/2 + lc - 1 <= nc
                R7 = corrcoef(A2(n+1:nt, ic-dc/2), A2(1:(nt-n), ic-dc/2+lc-1));
                R2(lc) = R7(2,1);
            else
                R2(lc) = NaN;
            end
        end
        R2M(ic - dc/2, :) = R2;
    end

    % Wavelength calculation
    L2M = zeros(1, size(R2M, 1));
    for i = 1:size(R2M, 1)
        try
            [~, ~, ~, trms, ~, ~, ~, ~] = Wave_Char(R2M(i,:), dt, 0, 2);
            L2M(i) = trms;
        catch ME
            fprintf('Error in Wave_Char (L2M) at i=%d: %s\n', i, ME.message);
            L2M(i) = NaN;
        end
    end

    % Wave period and Hs calculation
    T2M = nan(1, nc);
    Hs = nan(1, nc);
    for i = 1 + dc/2 : nc - dc/2
        try
            [hs, ~, ~, trms, ~, ~, ~, ~] = Wave_Char(detrend(A2(1:min([500 round(length(A2(:,1)))]), i)), dt, 0, 2);
            T2M(i) = trms;
            Hs(i) = hs;
        catch ME
            fprintf('Error in Wave_Char (T2M, Hs) at i=%d: %s\n', i, ME.message);
            T2M(i) = NaN;
            Hs(i) = NaN;
        end
    end

    RM = R2M;
    [~, R2M] = max(R2M');
    R2M = R2M ./ dpha;

    % Check for NaN values in outputs
    if any(isnan(R2M(:))) || any(isnan(L2M)) || any(isnan(T2M)) || any(isnan(Hs))
        warning('NaN values present in the output');
    end

catch ME
    % Error handling
    error('Error in CrossCorrelation_Coastcams: %s', ME.message);
end

end

