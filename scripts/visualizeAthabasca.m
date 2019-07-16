clear;

load('visdata.mat');

% Print it
clf; close all;
% Display scatter plots of gps surface data and bedrock picks
figure(1);
    subplot(2, 2, 1);
        scatter3(visdata.longs, visdata.lats, visdata.elevs, 10, visdata.elevs);
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Surface (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 2);
        scatter3(visdata.longPicks, visdata.latPicks, visdata.elevPicks-visdata.zPicks, 10, visdata.elevPicks-visdata.zPicks);
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Bedrock (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 3);
        scatter3(visdata.longsCSRS, visdata.latsCSRS, visdata.elevsCSRS, 10, visdata.elevsCSRS);
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Surface (CSRS data)");
        axis square;
        view(0, 90);
    subplot(2, 2, 4);
        scatter3(visdata.longCSRSPicks, visdata.latCSRSPicks, visdata.elevCSRSPicks-visdata.zPicks, 10, visdata.elevCSRSPicks-visdata.zPicks);
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("WIP Bedrock (CSRS data)");
        axis square;
        view(0, 90);
        
% Display surface mech plots of gps surface data and bedrock picks
figure(2);
    subplot(2, 2, 1);
        trisurf(visdata.tri, visdata.longs, visdata.lats, visdata.elevs)
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Surface (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 2);
        trisurf(visdata.triPick, visdata.longPicks, visdata.latPicks, visdata.elevPicks-visdata.zPicks)
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Bedrock (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 3);
        trisurf(visdata.triCSRS, visdata.longsCSRS, visdata.latsCSRS, visdata.elevsCSRS)
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("Surface (CSRS data)");
        axis square;
        view(0, 90);
    subplot(2, 2, 4);
        trisurf(visdata.triCSRSPick, visdata.longCSRSPicks, visdata.latCSRSPicks, visdata.elevCSRSPicks-visdata.zPicks)
        xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
        colorbar;
        title("WIP Bedrock (CSRS data)");
        axis square;
        view(0, 90);

% This figure displays downrange profiles of the glaciers
figure(3);
    % Sort it by the projection we'll display it in
    [~, sorted] = sort(visdata.longPicks-visdata.latPicks);
    [~, sortCSRS] = sort(visdata.longCSRSPicks-visdata.latCSRSPicks);
    % Select only those values on the downrange track we want
    index = abs(visdata.latPicks(sorted) + 0.777778*visdata.longPicks(sorted) - 143.389) < 0.0003;
    indexCSRS = abs(visdata.latCSRSPicks(sortCSRS) + 0.777778*visdata.longCSRSPicks(sortCSRS) - 143.389) < 0.0003;
    % Display it in map view with the picked values highlighted
    subplot(2, 2, 1); hold on;
        plot(visdata.latPicks, visdata.longPicks, '.');
        plot(visdata.latPicks(sorted(index)), visdata.longPicks(sorted(index)), 'o');
        title("Handheld GPS");
        axis square;
    subplot(2, 2, 2); hold on;
        plot(visdata.latCSRSPicks, visdata.longCSRSPicks, '.');
        plot(visdata.latCSRSPicks(sortCSRS(indexCSRS)), visdata.longCSRSPicks(sortCSRS(indexCSRS)), 'o');
        title("CSRS data");
        axis square;
    % Plot profiles
    subplot(2, 2, 3);
        plot(visdata.latPicks(sorted(index))-visdata.longPicks(sorted(index)), visdata.elevPicks(sorted(index))-visdata.zPicks(sorted(index)), 'o-');
        title("Handheld GPS");
    subplot(2, 2, 4);
        plot(visdata.latCSRSPicks(sortCSRS(indexCSRS))-visdata.longCSRSPicks(sortCSRS(indexCSRS)), visdata.elevCSRSPicks(sortCSRS(indexCSRS))-visdata.zPicks(sortCSRS(indexCSRS)), 'o-');
        title("CSRS data");        

% This figure serves as a legend for which gps data came from file
% figure(4); hold on;
%     N = length(bedrockPicks);
%     % Set up a data structure to build the legend off of
%     Legend = cell(N, 1);
%     % There's too many seperate traces to use color to tell them apart so
%     % we'll have to incorporate texture also
%     for i=1:N
%         switch mod(i, 4)
%             case 0
%                 linespec = '-';
%             case 1
%                 linespec = ':';
%             case 2
%                 linespec = '--';
%             case 3
%                 linespec = '-.';
%         end
%         d = bedrockPicks(i);
%         plot3(d.visdata.long, d.visdata.lat, d.visdata.elev-d.depth, linespec, 'LineWidth', 1.5, 'Color', hsv2rgb([floor(i/4)*4/N, 1, 0.8]));
%         plot3(d.visdata.long, d.visdata.lat, d.visdata.elev, linespec, 'LineWidth', 1.5, 'Color', hsv2rgb([floor(i/4)*4/N, 1, 0.8]));
%         Legend{i} = d.name;
%     end
%     % Add the highlighted values in the previous figure's profile
%     plot3(visdata.longPicks(sorted(index)), visdata.latPicks(sorted(index)), visdata.elevPicks(sorted(index))-visdata.zPicks(sorted(index)), 'o-', 'Color', [1, 0, 0]);
%     % Build legend
%     legend(Legend);
%     % Configure figure
%     xlabel("visdata.longitude"); ylabel("visdata.latitude"); zlabel("visdata.elevation (m)");
%     axis square;