%% Script to process Athabasca GPR data
%  06 May 2019
%  William Armstrong

clear all
close all

%% user inputs

dataFolderPath = 'C:\Users\armstrongwh\Google Drive\appstate\projects\athabasca\data\radar\martinRadar\res\';

% for radar gain enhancement
maxDepth = 400; % depth at which brightness factor will maximize [m]
maxBrighten = 7; % factor by which radar will be brightened at max depth
expDepthScale = 150; % length scale for exponential gain on radar gram
expPrefac = 1;

stdThresh = 2; % threshold for difference in standard deviation in difference of radargram in x direction (along profile) for defining stationary data


%% convert wv files to mat files and plot raw radargams

if 0 % turn on/off looping through wv folders to create mat files

    cd(dataFolderPath) % change path to data folder

    folders = dir; % get list of folders in diretory

    for i = 1:length(folders) % loop over folders in gpr data directory
        folderNow = folders(i); % get current folder

        if folderNow.isdir == 1 && folderNow.name(1) ~= '.' && strcmp(folderNow.name,'gps') ~= 1 % is file a folder, and not one of the root folders
            nameNow = folderNow.name;
            outputFn = [dataFolderPath,nameNow,'.mat']; % make name for output file

            if exist(outputFn) == 0 % if output file does not already exist
                rec = readwv(folderNow.name); % transform wv files to .mat file

                save(outputFn,'rec') % save output

            end

        end


    end

end

%% apply sine gain to radar data

cd(dataFolderPath) % change path to data folder

folders = dir([dataFolderPath '*.mat']); % get list of folders in diretory

for i = 1:length(folders)
    folderNow = folders(i);
    folderName = folderNow.name;
    
    try %  the strcmp errors with short filenames b/c not long enough for index, so use try to avoid this
        
        if strcmp(folderName(end-3:end),'.mat') ==1 && strcmp(folderName(end-7:end),'gain.mat') ~= 1 && strcmp(folderName(1:3),'mar') % this is gpr data file if true
            disp(['Processing: ',folderName])
            load(folderName)
            
            surveyName = folderName(1:end-4);
            
            recSine = rec; %copy rec file
            wvSize = size(recSine.wv); % dimensions of wv data
            depthLen = wvSize(2);
            shotLen = wvSize(1);
            
            % 30 m separation; 84 is the zero offset
            depth=(((1:1000)-84)*recSine.dt+30/300*1e-6)*169e6/2;
            
            maxDepthInd = max(find(depth<maxDepth)) ; % where is maximum depth for contrast enchangement located
            brightnessFactors = maxBrighten*sin( (depth/maxDepth) * (pi/2) ); % factor for enhancing contrast
            brightnessSquareFactors = maxBrighten*(depth/(maxDepth/2)).^2;
            brightnessExpFactors = expPrefac*exp(depth/(expDepthScale));
            
            recSine.wvSine = nan(shotLen,depthLen); % make empty storage matx
            recSine.wvSquare = nan(shotLen,depthLen); % make empty storage matx
            recSine.wvExp = nan(shotLen,depthLen); % make empty storage matx
            recSine.maxDepth = maxDepth; % store processing parameters so output file is stand-alone
            recSine.maxBrightenFactor = maxBrighten;
            
            % increase brightness at depth for each shot
            for j = 1:shotLen
                recSine.wvSine(j,:) = recSine.wv(j,:) .* brightnessFactors;
                recSine.wvSquare(j,:) = recSine.wv(j,:) .* brightnessSquareFactors;
                recSine.wvExp(j,:) = recSine.wv(j,:) .*brightnessExpFactors;             
            end
            
            % find stationary data
            
            dx = diff(recSine.wvSine,1,1); % compute difference in brightness in x direction (along profile)
            medDx = median(dx,2);
            stdDx = std(dx,0,2);
            
            movingIndDx = find(stdDx>stdThresh); % index for points where gpr is moving
            movingIndDx = [min(movingIndDx)-1; movingIndDx; max(movingIndDx)+1]; % add one buffer cell on either side
            movingInd = movingIndDx + 1; % for use on rec.wv (b/c difference makes one shorter)
            
            if 0
                figure(6)
                clf
                orient landscape
                plot(depth,brightnessExpFactors,'linewidth',4)
                xlim([0 400])
                set(gca,'fontsize',14,'linewidth',2)
                xlabel('Depth below surface [m]','fontsize',18)
                ylabel('Exponential gain correction factor [-]','fontsize',18)
                grid on
                
                h = gcf;
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);
                print('exponentialGainCorrection.pdf','-dpdf','-r300')
                print('exponentialGainCorrection.png','-dpng','-r300')                
                
                
                figure(5)
                clf
                orient landscape
                subplot(2,1,1)
                imagesc(1:shotLen-1,depth,recSine.wv')
                colormap('bone')
                caxis([-50,50]);
                ylim([0 350])

                subplot(2,1,2)
                imagesc(1:shotLen-1,depth,recSine.wvExp')
                colormap('bone')
                caxis([-50,50]);
                ylim([0 350])
            end
            
            if 0
                figure(4);
                clf
                orient landscape
                subplot(3,1,1:2)
                imagesc(1:shotLen-1,depth,dx');
                %colormap('bone')
                %colorbar()
                caxis([-50,50]);
                ylim([0 350])
                ylabel('Depth [m]','fontsize',18)
                %grid on
                set(gca,'fontsize',14,'linewidth',2)
                xl = xlim;
                yl = ylim;
                title([surveyName, ' data clipping; stdThresh = ', num2str(stdThresh)],'fontsize',24,'Interpreter', 'none')

                hold on
                plot([min(movingIndDx) min(movingIndDx)],yl,'k--','linewidth',2)
                plot([max(movingIndDx) max(movingIndDx)],yl,'k--','linewidth',2)

                subplot(3,1,3)
                plot(1:shotLen-1,stdDx,'linewidth',3);
                ylabel('\sigma_{dx}','fontsize',18)
                xlabel('Shot number [-]','fontsize',18)
                grid on
                set(gca,'fontsize',14,'linewidth',2)
                xlim(xl);
                yl = ylim;
                hold on
                plot([min(movingIndDx) min(movingIndDx)],yl,'k--','linewidth',2)
                plot([max(movingIndDx) max(movingIndDx)],yl,'k--','linewidth',2)
                plot(xl,[stdThresh,stdThresh],'k--','linewidth',2)

                h = gcf;
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);

                figOutFn = [surveyName,'_dataClipping.pdf'];
                pngOutFn = [surveyName,'_dataClipping.png'];

                %print(figOutFn,'-dpdf','-r300')
                %print(pngOutFn,'-dpng','-r300')
            end
            
            % save to recSine file
            recSine.movingInd = movingInd; 
            %recSine.wvSineClip = recSine.wvSine(movingInd,:);
            %recSine.wvClip = recSine.wv(movingInd,:);
            recSine.clipStdThresh = stdThresh;
            recSine.dxBrigthness = dx;
            recSine.stdDx = stdDx;
            recSine.wvExpDepthScale = expDepthScale;
            recSine.wvExpPrefactor = expPrefac;
            
            % save mat file out
            save([surveyName,'_gain.mat'],'recSine')
            
            if 1 % turns plotting and saving off/on
                
                % raw radar gram
                figure(1);
                clf
                orient landscape
                imagesc(1:shotLen,depth,recSine.wv');
                colormap('bone')
                caxis([-50,50]);
                ylim([0 350])
                
                yl = ylim;
                hold on
                plot([min(movingInd) min(movingInd)],yl,'r--','linewidth',2)
                plot([max(movingInd) max(movingInd)],yl,'r--','linewidth',2)
                
                set(gca,'fontsize',16,'linewidth',2)
                xlabel('Shot number [-]','fontsize',20)
                ylabel('Depth below surface [m]','fontsize',20)
                title([surveyName, ' raw radargram'],'fontsize',24,'Interpreter', 'none')

                h = gcf;
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);

                figOutFn = [surveyName,'_rawRadargram.pdf'];
                pngOutFn = [surveyName,'_rawRadargram.png'];

                %print(figOutFn,'-dpdf','-r300')
                %print(pngOutFn,'-dpng','-r300')



                % gain enhanced
                figure(2);
                clf
                orient landscape
                imagesc(1:shotLen,depth,recSine.wvExp');
                colormap('bone')
                caxis([-50,50]);
                ylim([0 350])
                yl = ylim;
                hold on
                plot([min(movingInd) min(movingInd)],yl,'r--','linewidth',2)
                plot([max(movingInd) max(movingInd)],yl,'r--','linewidth',2)                

                set(gca,'fontsize',16,'linewidth',2)
                xlabel('Shot number [-]','fontsize',20)
                ylabel('Depth below surface [m]','fontsize',20)
                title([surveyName, ' exp gain; depthScale = ',num2str(maxDepth) , ' m'],'fontsize',24,'Interpreter', 'none')

                h = gcf;
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);

                figOutFn = [surveyName,'_expGain.pdf'];
                pngOutFn = [surveyName,'_expGain.png'];

                print(figOutFn,'-dpdf','-r300')
                print(pngOutFn,'-dpng','-r300')


                % comparison figure
                figure(3);
                clf
                orient landscape
                subplot(2,1,1)
                imagesc(1:shotLen,depth,recSine.wv');
                colormap('bone')
                caxis([-50,50]);
                ylim([0 350])

                set(gca,'fontsize',16,'linewidth',2)
                xlabel('Shot number [-]','fontsize',20)
                ylabel('Depth below surface [m]','fontsize',20)
                title([surveyName ,' raw radargram'],'fontsize',20,'Interpreter', 'none')
                ylim([0 350])


                subplot(2,1,2)
                imagesc(1:shotLen,depth,recSine.wvExp');
                colormap('bone')
                %colorbar()
                caxis([-50,50]);
                ylim([0 350]) 

                set(gca,'fontsize',16,'linewidth',2)
                xlabel('Shot number [-]','fontsize',20)
                ylabel('Depth below surface [m]','fontsize',20)
                title(['exp gain; depthScale = ',num2str(expDepthScale) , ' m'],'fontsize',20,'Interpreter', 'none')
                ylim([0 350])

                h = gcf;
                set(h,'PaperUnits','normalized');
                set(h,'PaperPosition', [0 0 1 1]);

                figOutFn = [surveyName,'_rawGainEnhanceComparison.pdf'];
                pngOutFn = [surveyName,'_rawGainEnhanceComparison.png'];

                print(figOutFn,'-dpdf','-r300')
                print(pngOutFn,'-dpng','-r300')

            end
            
        else    
            disp(['Skipping: ',folderName])
        end
    
    catch % this will be run on short filenames (don't want to process these, they are not gpr data files)
        disp(['Skipping: ',folderName])
    end
    
end

%% plot x,y,z of aqcuisition

if 0 % turn on/off this block

    cd(dataFolderPath) % change path to data folder

    folders = dir([dataFolderPath,'*gain.mat']); % get list of *gain.mat folders in diretory

    % make empty storage matrices
    % 40 surveys, biggest survey has 560 shots
    timeMatx = nan(40,560);
    latMatx = nan(40,560);
    lonMatx = nan(40,560);
    elevMatx = nan(40,560);

    goodIterNum = 1;

    surveyNameList = [];

    for i = 1:length(folders)
        folderNow = folders(i);
        folderName = folderNow.name;

        disp([num2str(i), ' Processing: ',folderName])
        load(folderName)

        obsNum = length(recSine.time); % number of shots in survey

        timeMatx(goodIterNum,1:obsNum) = recSine.time; % append current time observations to list
        latMatx(goodIterNum,1:obsNum) = recSine.lat; % append current latitudes 
        lonMatx(goodIterNum,1:obsNum) = recSine.long;  % append current longitudes 
        elevMatx(goodIterNum,1:obsNum) = recSine.elev;  % append current longitudes

        % collect survey names
        surveyName = folderName(1:end-9);

        if length(surveyName) == 6 % to make all surveys same number of characters, needed for list of survey names
            surveyName = [surveyName(1:5) '0', surveyName(6)]; 
        end

        surveyNameList = [surveyNameList; surveyName];

        listLen = size(surveyNameList);


        goodIterNum = goodIterNum + 1; % update survey count

    end

end

%% flatten matrices to lists (i.e., size = 1,n)

if 0 % turn on/off this block

    matxSize = size(latMatx); % how many surveys are there?

    timeList = [];
    latList = [];
    lonList = [];
    elevList = [];
    profileNumList = [];
    newProfileInd = 1;
    uniqueNum = [];


    cumLen = 1; % for marking survey start locations 

    for i = 1:matxSize(1)
        timeList = [timeList, timeMatx(i,:)];
        latList = [latList, latMatx(i,:)];
        lonList = [lonList, lonMatx(i,:)];
        elevList = [elevList, elevMatx(i,:)];

        numRepeat = i*ones(1,length(timeMatx(i,:)));
        profileNumList = [profileNumList,numRepeat];

        cumLen = cumLen + length(timeMatx(i,:)); % how many items have been appended
        newProfileInd = [newProfileInd, cumLen];

        uniqueList = [newProfileInd(i):1:newProfileInd(i+1)-1];
        uniqueNum = [uniqueNum, uniqueList];

    end


    % save coordiantes out as a csv
    outData = [uniqueNum; profileNumList; timeList; latList; lonList; elevList];
    dlmwrite('radarCoordinatesAndTime.csv',outData','delimiter', ',', 'precision', 12)

end


%% plot lat lon data

if 0 % turn on/off this block

    matxSize = size(latMatx); % how many surveys are there?


    figure(1)
    clf
    orient landscape
    hold on


    for i = 1:matxSize(1) % iterate over all surveys

        lonNow = lonMatx(i,:);
        latNow = latMatx(i,:);

        goodDataInd = lonNow > 115 & lonNow < 120 & latNow > 51 & latNow < 54;
        goodDataNum = sum(goodDataInd);
        nanNum = sum(isnan(lonNow));
        badDataNum = length(lonNow) - nanNum - goodDataNum;

        goodLon = lonNow(goodDataInd);
        goodLat = latNow(goodDataInd);

        figure(1)
        scatter(goodLon,goodLat,20,'filled')
        %text(max(lonNow(goodDataInd)),max(latNow(goodDataInd)),surveyNameList(i,:))
        text(goodLon(end),goodLat(end),surveyNameList(i,:),'fontsize',12,'Interpreter', 'none')
        %pause(1)
    end


    set(gca,'fontsize',14,'linewidth',2,'xdir','reverse')
    xlabel('Longituide [^{\circ}W]','fontsize',18)
    ylabel('Latitude [^{\circ}N]','fontsize',18)
    axis equal
    grid on

    h = gcf;
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);

    figOutFn = ['surveys_latLon_withSurveyLabels.pdf'];
    pngOutFn = ['surveys_latLon_withSurveyLabels.png'];

    %print(figOutFn,'-dpdf','-r300')
    %print(pngOutFn,'-dpng','-r300')


end 
%% plot a single radargram

if 1 % turn on/off this block

    surveySineName = 'mar6_4_gain.mat';

    load(surveySineName)
    
    dataSize = size(recSine.wvSine);
    shotLen = dataSize(1);
    
    clipSize = size(recSine.wvClip);
    clipLen = clipSize(1);
    
    
    figure(1);
    clf
    orient landscape
    
    %imagesc(1:shotLen,depth,recSine.wvSine'); % gain enhanced no clip
    imagesc(1:clipLen,depth,recSine.wvSineClip'); % gain enhanced with clip
    colormap('bone')
    caxis([-50,50]);
    ylim([0 350])

    set(gca,'fontsize',16,'linewidth',2)
    xlabel('Shot number [-]','fontsize',20)
    ylabel('Depth below surface [m]','fontsize',20)
    %title(surveySineName,'fontsize',24,'Interpreter', 'none')

    h = gcf;
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    %print('exampleRadarData.png','-dpng','-r300')
    %print('exampleRadarData.pdf','-dpdf','-r300')

end
