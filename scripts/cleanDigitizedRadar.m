clear;

% Edit this variable to load other 
fileName = 'mar6_8';

% You should not need to change these
folderPath = '..\digitizedRadar\raw\';
outPath = '..\digitizedRadar\clean\';

% Load the data
load([folderPath, fileName, '_bedPicks.mat']);
d = recDigitize;

% These values configure the cleaning function used below
s = 5; % Sliding window size for gaussian convolution
b = 1; % Number of standard deviations to retain

% Make copies of all the data to be cleaned
d.latTrim = d.lat;
d.longTrim = d.long;
d.elevTrim = d.elev;
d.timeTrim = d.time;

% Impose additional bounds based on physical constraints
d.latTrim(d.lat<52 | d.lat>60) = NaN;
d.longTrim(d.long<100 | d.long>150) = NaN;
d.elevTrim(d.elev<1000 | d.elev>4000) = NaN;
d.timeTrim(d.time<7.35e5 | d.time>7.39e5) = NaN;

% Find outliers
d.latClean = removeOutliers(d.latTrim, b, s);
d.longClean = removeOutliers(d.longTrim, b, s);
d.elevClean = removeOutliers(d.elevTrim, b, s);
d.timeClean = removeOutliers(d.timeTrim, b, s);

% Get timespan limits
timemin = min(d.timeClean); timemax = max(d.timeClean);

% This is how many new values we're going to generate in our artifical time
% series
shotsno = length(d.timeClean); 

% Pretty straightforward, make shotsno number of evenly-spaced times from 
% the min through the max time
artificalTimeSeries = timemin:(timemax-timemin)/(shotsno-1):timemax; 

% This is an interpolation map, relating the old timevalues to their 
% location in the new timeseries.
d.interp = arrayfun(@(x) genInterpolationMap(x, d.timeClean), artificalTimeSeries);

% Interpolate everything using the interpolation map
d.latInterp = interpolateFromMap(d.interp, d.latClean);
d.longInterp = interpolateFromMap(d.interp, d.longClean);
d.elevInterp = interpolateFromMap(d.interp, d.elevClean);
d.timeInterp = interpolateFromMap(d.interp, d.timeClean);

% Get interpolated data about our picks from the values nearest to the
% picked times
d.latPicks = arrayfun(@(x) interpolatePick(x, d.latInterp), d.xPick); 
d.longPicks = arrayfun(@(x) interpolatePick(x, d.longInterp), d.xPick);
d.elevPicks = arrayfun(@(x) interpolatePick(x, d.elevInterp), d.xPick);
d.timePicks = arrayfun(@(x) interpolatePick(x, d.timeInterp), d.xPick);

% Display cleaning results
clf;
    % Latitude
    subplot(4,2,1); hold on; 
        % Plot values
        plot(d.timeTrim, d.latTrim, '*'); 
        plot(d.timeClean, d.latClean, '*'); 
        plot(d.timeInterp, d.latInterp, '.');
        % Recalculate the cutoff margins for display
        latSmooth = gaussSmooth(d.latTrim, s);
        latDiff = b * nanstd(d.latTrim-latSmooth);
        % Plot these cutoff margins
        plot(d.timeTrim, latSmooth, '--');
        plot(d.timeTrim, latSmooth - latDiff, '--'); 
        plot(d.timeTrim, latSmooth + latDiff, '--');
        % Configure
        ylim([min(d.latTrim), max(d.latTrim)]);
        title('Lat');
    hold off;

    % Longitude
    subplot(4,2,2); hold on; 
        % Plot values
        plot(d.timeTrim, d.longTrim, '*'); 
        plot(d.timeClean, d.longClean, '*'); 
        plot(d.timeInterp, d.longInterp, '.');
        % Recalculate the cutoff margins for display
        longSmooth = gaussSmooth(d.longTrim, s);
        longDiff = b * nanstd(d.longTrim-longSmooth);
        % Plot these cutoff margins
        plot(d.timeTrim, longSmooth, '--');
        plot(d.timeTrim, longSmooth - longDiff, '--'); 
        plot(d.timeTrim, longSmooth + longDiff, '--');
        % Configure
        ylim([min(d.longTrim), max(d.longTrim)]);
        title('Long');
    hold off;

    % Elevation
    subplot(4,2,3); hold on; 
        % Plot values
        plot(d.timeTrim, d.elevTrim, '*'); 
        plot(d.timeClean, d.elevClean, '*'); 
        plot(d.timeInterp, d.elevInterp, '.');
        % Recalculate the cutoff margins for display
        elevSmooth = gaussSmooth(d.elevTrim, s);
        elevDiff = b * nanstd(d.elevTrim-elevSmooth);
        % Plot these cutoff margins
        plot(d.timeTrim, elevSmooth, '--');
        plot(d.timeTrim, elevSmooth - elevDiff, '--'); 
        plot(d.timeTrim, elevSmooth + elevDiff, '--');
        % Configure
        ylim([min(d.elevTrim), max(d.elevTrim)]);
        title('Elev');
    hold off;

    % Time
    subplot(4,2,4); hold on; 
        % Plot values
        plot(d.timeTrim, d.timeTrim, '*'); 
        plot(d.timeClean, d.timeClean, '*'); 
        plot(d.timeInterp, d.timeInterp, '.');
        % Recalculate the cutoff margins for display
        timeSmooth = gaussSmooth(d.timeTrim, s);
        timeDiff = b * nanstd(d.timeTrim-timeSmooth);
        % Plot these cutoff margins
        plot(d.timeTrim, timeSmooth, '--');
        plot(d.timeTrim, timeSmooth - timeDiff, '--'); 
        plot(d.timeTrim, timeSmooth + timeDiff, '--');
        % Configure
        ylim([min(d.timeTrim), max(d.timeTrim)]);
        title('Time');
    hold off;

    % Visualize data in 3-space
    subplot(4,2,[5,6,7,8]); hold on;
        % Plot the surface 
        scatter3(d.latInterp, d.longInterp, d.elevInterp, 6, d.timeInterp, 'filled');
        % Plot the surface track of the bedrock picks
        scatter3(d.latPicks, d.longPicks, d.elevPicks, 30, d.timePicks); 
        % Plot the bedrock picks
        scatter3(d.latPicks, d.longPicks, d.elevPicks-d.zPick, 30, 'filled');
        % Configure
        xlabel('Lat'); ylabel('Long'); zlabel('Elev'); grid on; title('Map');
    hold off;

% Save it
recDigitize = d;
save([outPath fileName, '_cleanPicks.mat'], 'recDigitize')

% The smoothing function
function out = gaussSmooth (data, window)
    % Make a gaussian function
    kern = normpdf(1:window, 1+(window-1)/2, 1);
    kern = kern / sum(kern);
    % Set up the margins of the data to be convolved
    bottom = ones(1, floor(window/2)) * data(1);
    top = ones(1, floor(window/2)) * data(end);
    % Gaussian convolution
    out = conv([bottom, data, top], kern, 'valid');
end

% The cleaning function, b is the number of standard deviations to cutoff
% and s is the window with which to smooth the data
function out = removeOutliers (in, b, s)    
    % Setup data to output
    out = in;
    % Get smoothed version of the data
    M = gaussSmooth(in, s);
    % Calculate b number of standard deviations
    diff = b * nanstd(M-in);
    % Set outliers from said b number of standard deviations to NaNs
    out(in < M-diff) = NaN;
    out(in > M+diff) = NaN;
end

% Interpolation functions
function out = genInterpolationMap(val, data)
    j = 1;
    while isnan(data(j)) || data(j) < val
        j = j + 1;
    end
    up = j; upval = data(j);
    while isnan(data(j)) || data(j) > val
        j = j - 1;
    end
    down = j; downval = data(j);
    frac = 0;
    if upval-downval
        frac = (val-downval)/(upval-downval);
    end
    out = down + frac*(up-down);
end

function out = interpolateFromMap(interpMap, data)
    out = zeros(1, length(interpMap));
    for i = 1:length(interpMap)
        val = interpMap(i);
        j = 1;
        while isnan(data(j)) || j < val
            j = j + 1; 
            if j > length(data) 
                break; 
            end
        end
        if j > length(data) 
            out(i) = NaN; continue; 
        end
        up = j; upval = data(j);
        while isnan(data(j)) || j > val
            j = j - 1; 
            if ~j 
                break; 
            end
        end
        if ~j 
            out(i) = NaN; continue; 
        end
        down = j; downval = data(j);
        frac = 0;
        if up-down
            frac = (val-down)/(up-down);
        end
        out(i) = downval*(1-frac) + upval*frac;
    end
end


% The pick interpolation function
function out = interpolatePick (pick, data)
    d = [data(1), data, data(end)];
    int = floor(pick);
    frac = pick - int;
    out = d(int+1)*(1-frac) + d(int-1)*frac;
end