cd Data/AA_1__AB_20/Vasculature_Images
%imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = dir(fullfile('*.png'));
imageNames = {imageNames.name}';

L = dir;
A = {L.name};

dumbimage = {};

for ii = 10:10:300
    dumbimage = cat(2,dumbimage,strcat('PFig_',num2str(ii),'.png'));
end

outputVideo = VideoWriter(fullfile('video.mp4'),'MPEG-4');
outputVideo.FrameRate = 3;
open(outputVideo)

for ii = 1:length(dumbimage)
   img = imread(dumbimage{ii});
   writeVideo(outputVideo,img)
end

close(outputVideo)
cd ../../..