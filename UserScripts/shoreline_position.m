function [shoreline_positions] = shoreline_position(I, threshold, res, ShoreMethod, plotoption_shore)
% SHORELINE_POSITION Detect shoreline positions in a timestack image
%
% Inputs:
%   I - Timestack Image (RGB)
%   threshold - Threshold value for shoreline detection (used in method 1)
%   res - Resolution (meters per pixel)
%   ShoreMethod - Method for shoreline detection:
%       1 = Grayscale method (specifically designed for rocky platform)
%       2 = Red minus Blue channel method
%       3 = Color divergence method
%   plotoption_shore - Flag to plot results (1 = plot, 0 = don't plot)
%
% Output:
%   shoreline_positions - Detected shoreline positions

% Input validation
validateattributes(I, {'uint8', 'uint16', 'double'}, {'3d'}, 'shoreline_position', 'I');
validateattributes(threshold, {'numeric'}, {'scalar', 'positive'}, 'shoreline_position', 'threshold');
validateattributes(res, {'numeric'}, {'scalar', 'positive'}, 'shoreline_position', 'res');
validateattributes(ShoreMethod, {'numeric'}, {'scalar', 'integer', '>=', 1, '<=', 3}, 'shoreline_position', 'ShoreMethod');
validateattributes(plotoption_shore, {'numeric'}, {'scalar', 'binary'}, 'shoreline_position', 'plotoption_shore');

switch ShoreMethod
    case 1
        shoreline_positions = grayscale_method(I, threshold, plotoption_shore);
    case 2
        shoreline_positions = color_channel_divergence(I, plotoption_shore);
    case 3
        shoreline_positions = color_divergence_method(I, res, plotoption_shore);
    otherwise
        error('Invalid ShoreMethod specified');
end

end

function shoreline_positions = grayscale_method(I, threshold, plotoption_shore)
    gray = rgb2gray(I);
    gray_adj = imadjust(gray, stretchlim(gray), [0 1], 1);
    gray_av = mean(gray_adj, 1, 'omitnan');
    J = imrotate(gray_av, 90);
    
    [~, c] = size(J);
    xti = [1 c c 1];
    yti = [0, 0, ceil(size(J,1)), ceil(size(J,1))];
    
    P = improfile(J, xti, yti);
    F = flip(P);
    Ps = movmean(F, 20);
    m = 1:length(Ps);
    
    [mins, min_locs] = findpeaks(-Ps, 'MinPeakProminence', 8);
    [maxs, max_locs] = findpeaks(Ps, 'MinPeakProminence', 8);
    
    peaks = [0 0; maxs max_locs];
    troughs = [mins min_locs];
    
    sp = size(peaks, 1);
    st = size(troughs, 1);
    s = max(sp, st);
    dif = [[peaks; zeros(s-sp, 2)], [troughs; zeros(s-st, 2)]];
    difs = [dif sum(dif(:,[1,3]), 2)];
    
    idx = find(difs(:,5) > threshold, 1);
    if ~isempty(idx)
        shoreline_positions = difs(idx,4);
    else
        shoreline_positions = 0;
    end
    
    if plotoption_shore
        plot_shoreline(I, m, Ps, min_locs, max_locs);
    end
end

function shoreline_positions = color_channel_divergence(I, plotoption_shore)
    ccd_av = mean(I, 1, 'omitnan');
    J_ccd = imrotate(ccd_av, 90);
    
    RmB = double(J_ccd(:,:,1)) - double(J_ccd(:,:,3));
    [pdf_values, pdf_locs] = ksdensity(RmB(:));
    
    thresh_otsu = multithresh(RmB);
    I1 = find(pdf_locs < thresh_otsu);
    I2 = find(pdf_locs > thresh_otsu);
    [~, J1] = max(pdf_values(I1));
    [~, J2] = max(pdf_values(I2));
    
    thresh_weightings = [1/3 2/3];
    thresh = thresh_weightings(1)*pdf_locs(I1(J1)) + thresh_weightings(2)*pdf_locs(I2(J2));
    
    if plotoption_shore
        plot_color_channel_divergence(pdf_locs, pdf_values, I1, I2, J1, J2, thresh);
    end
    
    [C, ~] = contour(RmB, [thresh thresh]);
    II = find(C(1,:) == thresh);
    shoreline_positions = C(2, II);
end

function shoreline_positions = color_divergence_method(I, res, plotoption_shore)
    RMB_I = double(I(:,:,1)) - double(I(:,:,3));
    level = graythresh(RMB_I);
    bw_RMB = imbinarize(RMB_I, level);
    
    se = strel('disk', 5);
    bw_RMB = imopen(bw_RMB, se);
    
    shoreline = zeros(1, size(bw_RMB, 2));
    for col = 1:size(bw_RMB, 2)
        row = find(bw_RMB(:,col), 1, 'first');
        shoreline(col) = ifelse(isempty(row), NaN, row);
    end
    
    shoreline_x = (1:size(bw_RMB, 2)) * res;
    shoreline_positions = shoreline * res;
    
    if plotoption_shore
        plot_color_divergence(I, shoreline_x, shoreline_positions);
    end
end

function plot_shoreline(I, m, Ps, min_locs, max_locs)
    figure;
    imagesc(I);
    yyaxis right;
    hold on;
    plot(m, Ps, 'w');
    plot(min_locs, Ps(min_locs), 'rv', 'MarkerFaceColor', 'r');
    plot(max_locs, Ps(max_locs), 'rs', 'MarkerFaceColor', 'b');
    set(gca, 'fontsize', 14);
    title('Shoreline Detection - Grayscale Method');
    xlabel('Cross-shore distance (pixels)');
    ylabel('Intensity');
end

function plot_color_channel_divergence(pdf_locs, pdf_values, I1, I2, J1, J2, thresh)
    figure;
    plot(pdf_locs, pdf_values);
    hold on;
    plot(pdf_locs([I1(J1) I2(J2)]), pdf_values([I1(J1) I2(J2)]), 'ro');
    YL = ylim;
    plot([thresh thresh], YL, 'r:', 'linewidth', 2);
    xlabel('Red minus Blue', 'fontsize', 10);
    ylabel('Counts', 'fontsize', 10);
    title('Color Channel Divergence Method');
end

function plot_color_divergence(I, shoreline_x, shoreline_y)
    figure;
    imshow(I);
    hold on;
    plot(shoreline_x, shoreline_y, 'r', 'LineWidth', 2);
    title('Shoreline Detection - Color Divergence Method');
    xlabel('Cross-shore distance (m)');
    ylabel('Alongshore distance (m)');
end

function res = ifelse(condition, true_value, false_value)
    if condition
        res = true_value;
    else
        res = false_value;
    end
end