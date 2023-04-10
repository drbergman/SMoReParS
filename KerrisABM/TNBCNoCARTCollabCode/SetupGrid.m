%SetupGrid
%This sets up the initial Grid that consists of the voxels.  This is a
%structure of arrays with the different species and (possibly agents).
%Each row will contain one 'variable'.
%I got rid of Will's setup segment agents
%Changing how we set up the matrices - logical etc.

%Version 3
%Last revision: Kerri - April 16, 2012
%               Will - Feb 2, 2012
%               Kerri 3/13/13

%Params
initO2 = .02;
maxsegments = 10000;
gridsize = [500,500,500];
% for tumorsphere
VEGFmax = 30;
VEGFdenom = 1.67;

%Initialize the grid
voxelgrid = [];
VEGFgrad = {};

%Setup the Species
%voxelgrid.O2(1:gridsize(1),1:gridsize(2),1:gridsize(3)) = initO2;

%Setup the Agentgrid %ADDED BY KERRI
% voxelgrid.Agent(1:gridsize(1),1:gridsize(2),1:gridsize(3)) = 0;
% %Preallocate Indices
% voxelgrid.CapInd(1:gridsize(1),1:gridsize(2),1:gridsize(3)) = 0;
% voxelgrid.PosInd(1:gridsize(1),1:gridsize(2),1:gridsize(3)) = 0;

voxelgrid.Agent = false(gridsize(1),gridsize(2),gridsize(3));
voxelgrid.CapInd = zeros(gridsize(1),gridsize(2),gridsize(3), 'single');
voxelgrid.PosInd = zeros(gridsize(1),gridsize(2),gridsize(3),'single');

vector = zeros(5*5*5,3);
count = 1;
for x = -2:2
    for y=-2:2
        for z = -2:2
            vector(count,:) = [x y z];
            count = count+1;
        end
    end
end
vector2 = zeros(3*3*3,3);
count = 1;
for x=-1:1
    for y=-1:1
        for z=-1:1
            vector2(count,:) = [x y z];
            count = count+1;
        end
    end
end

voxelgrid.Check = setdiff(vector,vector2,'rows');  %This finds the values in vector not in vector2

ActiveList = zeros(maxsegments,1);

%% Vegf initialization CHOOSE MODE 1=Amina's 0=flat gradient w/ our test parameters
%VEGFmode = 'constant'; %'gradient','aminas', 'constant' %Kerri Now in
%Control

switch lower(VEGFmode)

    case 'aminas'
            for x=1:gridsize(1)
                for y=1:gridsize(2)
                    for z=1:gridsize(3)
                        C1 = 19;
                        C2 = 4;
                        C3 = .1;
                        default = .6;
                        if x>gridsize(1)/8 && x<gridsize(1)*.25 && z>gridsize(3)/4 && y>gridsize(2)*.25 && y<gridsize(2)*.75
                            mean = abs(x*C1/gridsize(3) - C2);
                            sigma = C3*mean;
                            voxelgrid.VEGF(x,y,z) = normrnd(mean, sigma); %mean; %mean for speed - normrnd(mean, sigma);
                        else 
                            voxelgrid.VEGF(x,y,z) = default;
                        end
                        
                    end
                end
            end
            
    case 'gradient'
                
            C1 = 1;
            [i,j,k] = meshgrid(1:gridsize(1),1:gridsize(2),1:gridsize(3)); %NOTE K: should depend on gridsize!!!
            [vx, vy, vz] = gradient(i.*C1/100 + .53);%gradient(((-(i-100).^2-(j-100).^2))./100 + .53);  %abs((i.^2-j.^2).*C1./100 + .53)
            size(vx)
            size(vy)
            size(vz)
            voxelgrid.VEGFgrad = zeros(gridsize(1),3,gridsize(2),gridsize(3)); %(x, gradient vector, y, z), this is because convenient to store vectors in the second dimension of the matrix

            for x=1:gridsize(1)
                for y=1:gridsize(2)
                    for z=1:gridsize(3)
                        C1 = 1;
                        C2 = .51;
                        C3 = .05;
                        mean = abs(x*C1/100 + .53); %(x.*y)./100 + .53; %abs(x*C1/100 + .53);
                        sigma = C3;
                        voxelgrid.VEGF(x,y,z) = mean; %mean for speed - normrnd(mean, sigma);
                        voxelgrid.VEGFgrad(x,:,y,z) = [vx(x,y,z) vy(x,y,z) vz(x,y,z)];
                    end
                end
            end
            
    case 'constant'
        uniform = 20;%.6; %Changed from .06 - Kerri ng/ml to match experiments
        voxelgrid.VEGFgrad = zeros(gridsize(1),3,gridsize(2),gridsize(3)); %(x, gradient vector, y, z), this is because convenient to store vectors in the second dimension of the matrix
        voxelgrid.VEGF = uniform*ones(gridsize(1),gridsize(2),gridsize(3));
        voxelgrid.VEGFgrad(:,:,:,:) = 0;
        
    case 'tumorsphere' %V2.5
        for x=1:gridsize(1)
                for y=1:gridsize(2)
                    for z=1:gridsize(3)
                        if ((x-50)^2 + (y-50)^2 + (z-50)^2 ) < 25^2
                            %disp('Its within the tumor')
                            voxelgrid.VEGF(x,y,z) = 15; %Highest conc
                        else
                            %disp('not in the tumor');
                            [theta phi r] = cart2sph(x-50,y-50,z-50);
                            tgrad = VEGFmax - (r/VEGFdenom);
                            %disp(tgrad)
                            voxelgrid.VEGF(x,y,z) = max(0, tgrad); %Highest conc
                        end
                    end
                end
        end
        
        %Make the gradient
        [px,py,pz] = gradient(voxelgrid.VEGF);
        for x = 1:gridsize(1)
            for y = 1:gridsize(2)
                for z = 1:gridsize(3)
                    voxelgrid.VEGFgrad(x,:,y,z) = [px(x,y,z) py(x,y,z) pz(x,y,z)];
                end
            end
        end
        figure
        %surf(px,py,pz)
        
        %To plot VEGF gradient
        [xm,ym,zm] = meshgrid(1:1:gridsize(1),1:1:gridsize(2),1:1:gridsize(3));
        slice(xm,ym,zm,voxelgrid.VEGF, [50 100], [50 100], [100 200]), shading flat

        
    end
        