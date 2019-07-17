clear;

% Load CSRS
csrsPath = '..\gps\';
csrsName = 'gps';
load([csrsPath, csrsName, '.mat']);

% Load other things
folderPath = '..\digitizedRadar\clean\';
files = dir([folderPath '*.mat']);

% Set up index variables
lats = []; longs = []; elevs = []; times = []; sources = [];
latPicks = []; longPicks = []; elevPicks = []; timePicks = []; sourcePicks = [];
latsCSRS = []; longsCSRS = []; elevsCSRS = []; timesCSRS = [];
zPicks = [];

sourceNames = {};

% Load base GPR stuff
i = 0;
for file = files'
    load([folderPath file.name]);
    i = i + 1;
    % Load the base GPR data
    lats = [lats, recDigitize.latInterp];
    longs = [longs, recDigitize.longInterp];
    elevs = [elevs, recDigitize.elevInterp];
    times = [times, recDigitize.timeInterp];
    sources = [sources, i*ones(1, length(recDigitize.timeInterp))];
    % Load the shot pick data
    latPicks = [latPicks, recDigitize.latPicks'];
    longPicks = [longPicks, recDigitize.longPicks'];
    elevPicks = [elevPicks, recDigitize.elevPicks'];
    timePicks = [timePicks, recDigitize.timePicks'];
    sourcePicks = [sourcePicks, i*ones(1, length(recDigitize.timePicks))];
    % Save depths
    zPicks = [zPicks, recDigitize.zPick'];     
    % Save source names
    sourceNames{i} = file.name;
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
sources(nans == 1) = [];
latPicks(nanPicks == 1) = [];
longPicks(nanPicks == 1) = [];
elevPicks(nanPicks == 1) = [];
timePicks(nanPicks == 1) = [];
sourcePicks(nanPicks == 1) = [];
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
sources = sources(index);
latPicks = latPicks(pickindex);
longPicks = longPicks(pickindex);
elevPicks = elevPicks(pickindex);
timePicks = timePicks(pickindex);
sourcePicks = sourcePicks(pickindex);
zPicks = zPicks(pickindex);
latsCSRS = latsCSRS(CSRSindex);
longsCSRS = longsCSRS(CSRSindex);
elevsCSRS = elevsCSRS(CSRSindex);
timesCSRS = timesCSRS(CSRSindex);

% Get time indices for CSRS correlation
timeCSRSPicks = timePicks;

% Copy over sources
sourceCSRSPicks = sourcePicks;

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
visdata = {};
visdata.times = times;
visdata.lats = lats;
visdata.longs = longs;
visdata.elevs = elevs;
visdata.tri = tri;
visdata.sources = sources;
visdata.timePicks = timePicks;
visdata.latPicks = latPicks;
visdata.longPicks = longPicks;
visdata.elevPicks = elevPicks;
visdata.triPick = triPick;
visdata.sourcePicks = sourcePicks;
visdata.timesCSRS = timesCSRS;
visdata.latsCSRS = latsCSRS;
visdata.longsCSRS = longsCSRS;
visdata.elevsCSRS = elevsCSRS;
visdata.triCSRS = triCSRS;
visdata.timeCSRSPicks = timeCSRSPicks;
visdata.latCSRSPicks = latCSRSPicks;
visdata.longCSRSPicks = longCSRSPicks;
visdata.elevCSRSPicks = elevCSRSPicks;
visdata.triCSRSPick = triCSRSPick;
visdata.sourceCSRSPicks = sourceCSRSPicks;
visdata.zPicks = zPicks;
visdata.sourceNames = sourceNames;

save('visdata.mat', 'visdata')