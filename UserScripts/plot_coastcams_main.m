function plot_coastcams_main(Time_TS, Stack_av, SLA_S, Hs_TS, Tp_TS, rotation)
    % Create a time vector for SLA_S
    time_SLA = linspace(Time_TS(1), Time_TS(end), size(SLA_S, 1));

    % Create figure
    figure('Position', [100, 100, 1200, 800])

    % Subplot 1: Average Timestack
    ax1 = subplot(4,1,1);
    rotated_stack = imrotate(uint8(Stack_av), rotation);
    imagesc(Time_TS, 1:size(rotated_stack, 1), rotated_stack)
    colormap(ax1, 'gray')
    title('Average Timestack', 'FontSize', 14)
    ylabel('Cross-shore distance', 'FontSize', 12)
    set(ax1, 'YDir', 'normal')

    % Subplot 2: Sea Level Anomaly (SLA)
    ax2 = subplot(4,1,2);
    imagesc(time_SLA, 1:size(SLA_S, 2), SLA_S')
    colorbar
    title('Sea Level Anomaly (SLA)', 'FontSize', 14)
    ylabel('Timestack Image Number', 'FontSize', 12)
    colormap(ax2, 'jet')
    set(ax2, 'YDir', 'normal')

    % Subplot 3: Significant Wave Height
    ax3 = subplot(4,1,3);
    plot(Time_TS, Hs_TS, 'r.-')
    title('Significant Wave Height', 'FontSize', 14)
    ylabel('H_s [m]', 'FontSize', 12)
    grid on

    % Subplot 4: Peak Wave Period
    ax4 = subplot(4,1,4);
    plot(Time_TS, Tp_TS, 'b.-')
    title('Peak Wave Period', 'FontSize', 14)
    ylabel('T_p [s]', 'FontSize', 12)
    xlabel('Time', 'FontSize', 12)
    grid on

    % Link x-axes of all subplots
    linkaxes([ax1, ax2, ax3, ax4], 'x')

    % Adjust time axis labels
    datetick(ax1, 'x', 'HH:MM', 'keepticks', 'keeplimits')
    datetick(ax2, 'x', 'HH:MM', 'keepticks', 'keeplimits')
    datetick(ax3, 'x', 'HH:MM', 'keepticks', 'keeplimits')
    datetick(ax4, 'x', 'HH:MM', 'keepticks', 'keeplimits')
    
    % Adjust spacing between subplots
    set(gcf, 'Color', 'w')
end