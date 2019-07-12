%% Script for digitizing Athabasca bed returns
%  03 June 2019
%  W Armstrong

clear all
close all

%% load data

folderPath = 'C:\Users\hengstam\Documents\athabasca\radar\mat_files\gain\'; % folder path where *_gain.mat files live
outputPath  = 'C:\Users\hengstam\Documents\athabasca\digitizedRadar\';
profileName = 'mar7_10'; % profile to digitize

load([folderPath,profileName,'_gain.mat'])


%% plot

recDigitize = recSine; % copy file structure
recDigitize.profileName = profileName;

wvSize = size(recDigitize.wv); % dimensions of wv data
depthLen = wvSize(2);
shotLen = wvSize(1);

% 30 m separation; 84 is the zero offset
depth=(((1:1000)-84)*recDigitize.dt+30/300*1e-6)*169e6/2;

figure(1);
clf
orient landscape
set(gcf, 'Position', get(0, 'Screensize'));
imagesc(1:shotLen,depth,recDigitize.wvExp');
colormap('bone')
caxis([-50,50]);
ylim([0 350])
yl = ylim;
hold on            

set(gca,'fontsize',16,'linewidth',2)
xlabel('Shot number [-]','fontsize',20)
ylabel('Depth below surface [m]','fontsize',20)
title(profileName,'fontsize',24,'Interpreter', 'none')

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

[x,z] = ginput(1e3); % record cursor clicks; press enter to kill

figure(1)
hold on
scatter(x,z,100,'rx','linewidth',2)

% save digitized bed lcoations
recDigitize.xPick = x;
recDigitize.zPick = z;

% save figure out
figOutFn = [outputPath, profileName,'_bedPicks.pdf'];
pngOutFn = [outputPath, profileName,'_bedPicks.png'];
print(figOutFn,'-dpdf','-r300')
print(pngOutFn,'-dpng','-r300')

% save mat file out
save([outputPath, profileName, '_bedPicks.mat'],'recDigitize')

