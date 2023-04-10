%DoCoolPlot2Stem
%Just plot the stem cells
%1/31/13 KAN

iterstart = 100;
iterstep = 10;
iterend = 300;

AA = 75;
AB = 150;
BB = 2;
gradated=0;

%[tumor cells, vasc, CART, macrophages, hypoxic cells, antigens?]
params = [0,1,0,0,0,1];

movname = strcat('AA_',num2str(AA),'__AB_',num2str(AB),'BB_',num2str(BB),'_gradated',num2str(gradated));

spin = false;%true;
spinperiterstep = 360;
spinstep = 3;


if spin
    allTheFrames = cell(((iterend-iterstart)/iterstep)*(spinperiterstep/spinstep),1);
    allTheColorMaps = cell(((iterend-iterstart)/iterstep)*(spinperiterstep/spinstep),1);
    allTheColorMaps(:) = {zeros(256, 3)};
else
    allTheFrames = cell((iterend-iterstart)/iterstep,1);
    allTheColorMaps = cell((iterend-iterstart)/iterstep,1);
    allTheColorMaps(:) = {zeros(256, 3)};
end

myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
i = 1;

for iter = iterstart:iterstep:iterend
    
   f = PlotEverything(iter,AA,AB,BB,gradated,params(1),params(2),params(3),params(4),params(5),params(6));
   set(f,'Position',[10,10,1000,1000]);
   if ~spin
    myMovie(i) = getframe(f);
    i=i+1;
    %frame = getframe(f);
    %myMovie(iter).cdata = frame.cdata;
    %myMovie(iter).colormap = frame.colormap;
    
    %Save figure
    filenamep = ['seed',num2str(seedrate2), 'sym', num2str(100*symchange),'Mig', num2str(DoMigration),'Hyp', '_ii',num2str(ii), '.jpeg'];
    %eval(['print -dpict ',filenamep]); % following uses low bit graphic: imwrite(xx,mm,filename)
    print( '-dbmp', filenamep); %I think he's using a mac file
   else
       for j = 1:spinstep:spinperiterstep
           view(j,35);
           myMovie(i) = getframe(f);
           i=i+1;
       end
   end
   close(f);
end
cd Data
writerObj = VideoWriter(movname,'MPEG-4');
writerObj.Quality = 100;
if spin
   writerObj.FrameRate = 24;
else
writerObj.FrameRate = 5;
end
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(myMovie)
    % convert the image to a frame 
    writeVideo(writerObj,  myMovie(i));
end
% close the writer object
close(writerObj);
cd ..
