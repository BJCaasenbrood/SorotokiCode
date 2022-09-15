function gif2png(file,N)
% Pedro.F Rodenas Perez
% Get images from gif and save then to the folder

infos = imfinfo(file, 'gif');

[allframedata, ~] = imread(file, 'frames', 'all');
alldimensions = size(allframedata);
N0 = alldimensions(end);

setImg = round(linspace(1,N0,N));
j = 0;
for i=setImg
    [im, map] = imread(file, 'frames',i);
    %im = allframedata(:,:,1,i);
    num = num2str(j);
    % Tell your image names in first ''
    M = strcat('im',num,'.jpg');
    imwrite(im,map,M);
    
    j = j +1;
end

end

