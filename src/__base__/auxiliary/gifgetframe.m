function im = gifgetframe(file,N)
% Pedro.F Rodenas Perez
% Get images from gif and save then to the folder

%infos = imfinfo(file, 'gif');

%[allframedata, ~] = imread(file, 'frames', 'all');
%alldimensions = size(allframedata);
%N0 = alldimensions(end);

%if N < N0
    [im, map] = imread(file, 'frames',N);
%else
%    [im, map] = imread(file, 'frames',N);
%end

end

