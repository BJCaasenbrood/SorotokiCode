function jpg2mov(workingDir,file,fps)
%JPG2MOV converts a folder of jpg files into mp4-formatted video
imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';

tmp = cellfun(@(x) erase(x,file),imageNames,'UniformOutput',false);
numNames = cellfun(@(x) str2double(erase(x,'.png')),tmp,'UniformOutput',false);

[~,I] = sort(vertcat(numNames{:}));
outputVideo = VideoWriter(fullfile(workingDir,'tmp.mp4'));
outputVideo.FrameRate = fps;
open(outputVideo)

for ii = 1:length(I)
   [x,map] = imread(fullfile(workingDir,imageNames{I(ii)}),'png');
   img = ind2rgb(x, map);

   writeVideo(outputVideo,img);
   cout('written frame: ');
   cout([num2str(ii),'\n']);
end

end

