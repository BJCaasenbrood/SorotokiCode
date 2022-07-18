function gif2png(file,N)

% Pedro.F Rodenas Perez
% Get images from gif and save then to the folder

infos = imfinfo(file, 'gif');
%[im, map] = imread([GIF_file.folder, '\', GIF_file.name], 'frames', num);
%imshow(im)
% [allframedata, map] = imread(file, 'frames', 'all');
% alldimensions = size(allframedata);
% N = alldimensions(end);

for i=1:N
    [im, map] = imread(file, 'frames',i);
    %im = allframedata(:,:,1,i);
    num = num2str(i);
    % Tell your image names in first ''
    M = strcat('im',num,'.jpg');
    imwrite(im,map,M);
end
end

