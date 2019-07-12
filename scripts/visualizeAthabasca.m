clear;

% Load CSRS
csrsPath = 'C:\Users\hengstam\Documents\athabasca\gps\';
csrsName = 'gps';
load([csrsPath, csrsName, '.mat']);

% Load other things
files = dir('*.mat');

% Set up index variables
lats = []; longs = []; elevs = []; times = [];
latPicks = []; longPicks = []; elevPicks = []; timePicks = []; zPicks = [];
latsCSRS = []; longsCSRS = []; elevsCSRS = []; timesCSRS = [];

% Set up another data structure, used for the bedrock trace legend
bedrockPicks = struct(...
    'name', {},...
    'lat', {},...
    'long', {},...
    'elev', {},...
    'depth', {},...
    'time', {}...
);

% Load base GPR stuff
for file = files'
    load(file.name);
    % Load the base GPR data
    lats = [lats, recDigitize.latInterp];
    longs = [longs, recDigitize.longInterp];
    elevs = [elevs, recDigitize.elevInterp];
    times = [times, recDigitize.timeInterp];
    % Load the shot pick data
    latPicks = [latPicks, recDigitize.latPicks'];
    longPicks = [longPicks, recDigitize.longPicks'];
    elevPicks = [elevPicks, recDigitize.elevPicks'];
    timePicks = [timePicks, recDigitize.timePicks'];
    zPicks = [zPicks, recDigitize.zPick'];
    % Fill up the ground trace structure
    bedrockPicks(end+1) = struct(...
        'name', file.name(1:end-15),...
        'lat', recDigitize.latPicks',...
        'long', recDigitize.longPicks',...
        'elev', recDigitize.elevPicks',...
        'depth', recDigitize.zPick',...
        'time', recDigitize.timePicks');
        
end

% Load all CSRS data
for i=1:length(gps)
    % Load the CSRS data
    latsCSRS = [latsCSRS; gps{1, i}.latitude_decimal_degree];
    longsCSRS = [longsCSRS; -gps{1, i}.longitude_decimal_degree];
    elevsCSRS = [elevsCSRS; gps{1, i}.ellipsoidal_height_m];
    timesCSRS = [timesCSRS; datenum(2019, 3, gps{1, i}.day_of_year-59)+(gps{1, i}.decimal_hour/24)];
end

% NaN chart
nans = isnan(lats) | isnan(longs) | isnan(elevs) | isnan(times);
nanPicks = isnan(latPicks) | isnan(longPicks) | isnan(elevPicks) | isnan(timePicks);

% Remove NaNs
lats(nans == 1) = [];
longs(nans == 1) = [];
elevs(nans == 1) = [];
times(nans == 1) = [];
latPicks(nanPicks == 1) = [];
longPicks(nanPicks == 1) = [];
elevPicks(nanPicks == 1) = [];
timePicks(nanPicks == 1) = [];
zPicks(nanPicks == 1) = [];

% Get sorting index
[~, index] = sort(times);
[~, pickindex] = sort(timePicks);
[~, CSRSindex] = sort(timesCSRS);

% Sort everything
lats = lats(index);
longs = longs(index);
elevs = elevs(index);
times = times(index);
latPicks = latPicks(pickindex);
longPicks = longPicks(pickindex);
elevPicks = elevPicks(pickindex);
timePicks = timePicks(pickindex);
zPicks = zPicks(pickindex);
latsCSRS = latsCSRS(CSRSindex);
longsCSRS = longsCSRS(CSRSindex);
elevsCSRS = elevsCSRS(CSRSindex);
timesCSRS = timesCSRS(CSRSindex);

% Get time indices for CSRS correlation
timeCSRSPicks = timePicks;

% Repick positions based on CSRS stuff
latCSRSPicks = ones(1, length(timePicks));
longCSRSPicks = ones(1, length(timePicks));
elevCSRSPicks = ones(1, length(timePicks));
interp = ones(1, length(timePicks));
for i=1:length(timePicks) 
    interp(i) = 0;
    t = timeCSRSPicks(i);
    j = 1;
    while timesCSRS(j) < t
        j = j+1;
    end
    if j > 1
        if timesCSRS(j) < 1
            interp(i) = (t-timesCSRS(j-1))/(timesCSRS(j)-timesCSRS(j-1));
        end
        latCSRSPicks(i) = latsCSRS(j-1)*(1-interp(i)) + latsCSRS(j)*interp(i);
        longCSRSPicks(i) = longsCSRS(j-1)*(1-interp(i)) + longsCSRS(j)*interp(i);
        elevCSRSPicks(i) = elevsCSRS(j-1)*(1-interp(i)) + elevsCSRS(j)*interp(i);
    else 
        latCSRSPicks(i) = latsCSRS(1);
        longCSRSPicks(i) = longsCSRS(1);
        elevCSRSPicks(i) = elevsCSRS(1);
    end        
    interp(i) = interp(i)+j-1;
end

% Visualize data in 3-space using a DeLaunay triangular mesh
tri = delaunay(longs, lats);
triPick = delaunay(longPicks, latPicks);
triCSRS = delaunay(longsCSRS, latsCSRS);
triCSRSPick = delaunay(longCSRSPicks, latCSRSPicks);

% Save it all
% visdata = {};
% visdata.times = times;
% visdata.lats = lats;
% visdata.longs = longs;
% visdata.elevs = elevs;
% visdata.tri = tri;
% visdata.timePicks = timePicks;
% visdata.latPicks = latPicks;
% visdata.longPicks = longPicks;
% visdata.elevPicks = elevPicks;
% visdata.triPick = triPick;
% visdata.timesCSRS = timesCSRS;
% visdata.latsCSRS = latsCSRS;
% visdata.longsCSRS = longsCSRS;
% visdata.elevsCSRS = elevsCSRS;
% visdata.triCSRS = triCSRS;
% visdata.timeCSRSPicks = timeCSRSPicks;
% visdata.latCSRSPicks = latCSRSPicks;
% visdata.longCSRSPicks = longCSRSPicks;
% visdata.elevCSRSPicks = elevCSRSPicks;
% visdata.triCSRSPick = triCSRSPick;
% visdata.zPicks = zPicks;

% Print it
clf; close all;
% Display scatter plots of gps surface data and bedrock picks
figure(1);
    subplot(2, 2, 1);
        scatter3(longs, lats, elevs, 10, elevs);
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Surface (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 2);
        scatter3(longPicks, latPicks, elevPicks-zPicks, 10, elevPicks-zPicks);
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Bedrock (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 3);
        scatter3(longsCSRS, latsCSRS, elevsCSRS, 10, elevsCSRS);
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Surface (CSRS data)");
        axis square;
        view(0, 90);
    subplot(2, 2, 4);
        scatter3(longCSRSPicks, latCSRSPicks, elevCSRSPicks-zPicks, 10, elevCSRSPicks-zPicks);
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("WIP Bedrock (CSRS data)");
        axis square;
        view(0, 90);
        
% Display surface mech plots of gps surface data and bedrock picks
figure(2);
    subplot(2, 2, 1);
        trisurf(tri, longs, lats, elevs)
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Surface (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 2);
        trisurf(triPick, longPicks, latPicks, elevPicks-zPicks)
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Bedrock (Handheld GPS)");
        axis square;
        view(0, 90);
    subplot(2, 2, 3);
        trisurf(triCSRS, longsCSRS, latsCSRS, elevsCSRS)
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("Surface (CSRS data)");
        axis square;
        view(0, 90);
    subplot(2, 2, 4);
        trisurf(triCSRSPick, longCSRSPicks, latCSRSPicks, elevCSRSPicks-zPicks)
        xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
        colorbar;
        title("WIP Bedrock (CSRS data)");
        axis square;
        view(0, 90);

% This figure displays downrange profiles of the glaciers
figure(3);
    % Sort it by the projection we'll display it in
    [~, sorted] = sort(longPicks-latPicks);
    [~, sortCSRS] = sort(longCSRSPicks-latCSRSPicks);
    % Select only those values on the downrange track we want
    index = abs(latPicks(sorted) + 0.777778*longPicks(sorted) - 143.389) < 0.0003;
    indexCSRS = abs(latCSRSPicks(sortCSRS) + 0.777778*longCSRSPicks(sortCSRS) - 143.389) < 0.0003;
    % Display it in map view with the picked values highlighted
    subplot(2, 2, 1); hold on;
        plot(latPicks, longPicks, '.');
        plot(latPicks(sorted(index)), longPicks(sorted(index)), 'o');
        title("Handheld GPS");
        axis square;
    subplot(2, 2, 2); hold on;
        plot(latCSRSPicks, longCSRSPicks, '.');
        plot(latCSRSPicks(sortCSRS(indexCSRS)), longCSRSPicks(sortCSRS(indexCSRS)), 'o');
        title("CSRS data");
        axis square;
    % Plot profiles
    subplot(2, 2, 3);
        plot(latPicks(sorted(index))-longPicks(sorted(index)), elevPicks(sorted(index))-zPicks(sorted(index)), 'o-');
        title("Handheld GPS");
    subplot(2, 2, 4);
        plot(latCSRSPicks(sortCSRS(indexCSRS))-longCSRSPicks(sortCSRS(indexCSRS)), elevCSRSPicks(sortCSRS(indexCSRS))-zPicks(sortCSRS(indexCSRS)), 'o-');
        title("CSRS data");        

% This figure serves as a legend for which gps data came from file
figure(4); hold on;
    N = length(bedrockPicks);
    % Set up a data structure to build the legend off of
    Legend = cell(N, 1);
    % There's too many seperate traces to use color to tell them apart so
    % we'll have to incorporate texture also
    for i=1:N
        switch mod(i, 4)
            case 0
                linespec = '-';
            case 1
                linespec = ':';
            case 2
                linespec = '--';
            case 3
                linespec = '-.';
        end
        d = bedrockPicks(i);
        plot3(d.long, d.lat, d.elev-d.depth, linespec, 'LineWidth', 1.5, 'Color', hsv2rgb([floor(i/4)*4/N, 1, 0.8]));
        plot3(d.long, d.lat, d.elev, linespec, 'LineWidth', 1.5, 'Color', hsv2rgb([floor(i/4)*4/N, 1, 0.8]));
        Legend{i} = d.name;
    end
    % Add the highlighted values in the previous figure's profile
    plot3(longPicks(sorted(index)), latPicks(sorted(index)), elevPicks(sorted(index))-zPicks(sorted(index)), 'o-', 'Color', [1, 0, 0]);
    % Build legend
    legend(Legend);
    % Configure figure
    xlabel("Longitude"); ylabel("Latitude"); zlabel("Elevation (m)");
    axis square;