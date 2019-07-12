function rec=readwv(folder)

% function to read radar wave files and save as mat file
% input is folder name

files=dir([folder,'/*.txt']);
cd(folder);

for fn=1:size(files,1)

    fid = fopen(['wv' num2str(fn) '.txt']);
    nline=fgetl(fid);
    nnline=fgetl(fid);
    timestr=[nline, ' ', nnline];
    rec.time(fn)=datenum(timestr);

    % lat long elev

    latstr=fgetl(fid);
    latnum=str2num(latstr);
    latdeg=floor(latnum/100);
    rec.lat(fn)=latdeg+(latnum-latdeg*100)/60;

    longstr=fgetl(fid);
    longnum=str2num(longstr);
    longdeg=floor(longnum/100);
    rec.long(fn)=longdeg+(longnum-longdeg*100)/60;

    rec.elev(fn) = str2num(fgetl(fid));

    rec.dt = str2num(fgetl(fid));

    % read waveform
    nline=fgetl(fid);
    rec.wv(fn,:)=str2num(nline);
    fclose(fid);
    
end

% 30 m separation; 84 is the zero offset

depth=(((1:1000)-84)*rec.dt+30/300*1e-6)*169e6/2;

figure(1);
clf
orient landscape
imagesc(1:fn,depth,rec.wv');
colormap('bone')
caxis([-50,50]);
set(gca,'fontsize',16,'linewidth',2)
xlabel('Shot number [-]','fontsize',20)
ylabel('Depth below surface [m]','fontsize',20)
title(folder,'fontsize',24)
ylim([0 350])
cd ..

h = gcf;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);

figOutFn = [folder,'_rawRadargram.pdf'];
pngOutFn = [folder,'_rawRadargram.png'];

print(figOutFn,'-dpdf','-r300')
print(pngOutFn,'-dpng','-r300')

