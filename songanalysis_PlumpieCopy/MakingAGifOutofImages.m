function [] = MakingAGifOutofImages(ImageFileDir, ImageRoot, StartingFrameNum, EndingFrameNum)

% First change to the directory and then read in one file at a time and
% write to an image

PresentDir = pwd;
cd(ImageFileDir);
ImageFiles = dir([ImageRoot, '*']);
for i = 1:length(ImageFiles),
    ImageFileNames{i} = ImageFiles(i).name;
end

FrameStep = 3;
for i = StartingFrameNum:3:EndingFrameNum,
    FileIndex = find(cellfun(@length, strfind(ImageFileNames, num2str(i))));
    if (~isempty(FileIndex))
        FileIndex = FileIndex(1);
        [x] = imread(ImageFileNames{FileIndex});
        x = imresize(x, 0.3);
        [X, map] = rgb2ind(x, 16);
        if (i == StartingFrameNum)
            imwrite(X, map, [ImageRoot, '.Frame.', num2str(StartingFrameNum), '.to.', num2str(EndingFrameNum), '.gif'], 'Loopcount', inf, 'DelayTime', (FrameStep+1)/25);
        else
            imwrite(X, map, [ImageRoot, '.Frame.', num2str(StartingFrameNum), '.to.', num2str(EndingFrameNum), '.gif'], 'WriteMode', 'append', 'DelayTime', (FrameStep+1)/25);
        end
    end
end

disp('Finished');