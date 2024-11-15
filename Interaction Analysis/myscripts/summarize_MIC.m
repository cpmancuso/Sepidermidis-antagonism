function fighandle = summarize_MIC(filepath,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 700 400];

range = 'C1:L16'; 
[data, headers] = xlsread(filepath, range);

isolates = data(:,1);
data = data(:,2:end);
headers = headers(2:end);

numColors = 7;
cmap = colormap([linspace(1,1,numColors); linspace(0,1,numColors); linspace(0,1,numColors)].');

% Create heatmaps
for f=1:numel(headers)
    array = data(:,f);
    subplot(1,numel(headers),f)
    imagesc(log2(array));
    for n=1:numel(array)
        text(1,n,num2str(array(n)),'HorizontalAlignment','center')
    end
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title(headers{f})
end