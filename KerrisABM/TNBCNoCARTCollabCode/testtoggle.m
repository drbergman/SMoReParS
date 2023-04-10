function varargout = testtoggle(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @testtoggle_OpeningFcn, ...
                   'gui_OutputFcn',  @testtoggle_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before testtoggle is made visible.
function testtoggle_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to testtoggle (see VARARGIN)

% Choose default command line output for testtoggle
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes testtoggle wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = testtoggle_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value=get(hObject, 'Value');

if(value==1)
iter = 6;

%cd Data/AA_0__BB_5/HypLocations
cd Data/AA_300__BB_2/HypLocations
XYZ = dlmread('XYZloc_t_10.txt');
XYZkt = dlmread('XYZkt_t_10.txt');
seedrate2 = 0;
symchange = 0.02;
DoMigration = 3;
ii = 1080;
load('voxelgrid_t_10.mat');

plane = find(XYZ(:,3) == 25);

%XYZ = gpuArray(XYZHyp(:,:));
 %Stemh = find(Stateall(:,1) == 1);
 %CCR5h = find(Stateall(:,3) == 1);
 %both = find(Stateall(:,1) == 1 & Stateall(:,3) == 1);

 %Stateall(CCR5h,1) = 3;
%Stateall(both,1) = 4;


% xmin = 0; % need better estimates for these
% xmax = 500;  
% ymin = 0;  
% ymax = 500; 
firstState = 0;
numStates = 3;
         
voxelgrid.Agent(2,:,:) = voxelgrid.Agent(1,:,:)|voxelgrid.Agent(2,:,:);
voxelgrid.Agent(:,2,:) = voxelgrid.Agent(:,1,:)|voxelgrid.Agent(:,1,:);         
voxelgrid.Agent(:,:,2) = voxelgrid.Agent(:,:,1)|voxelgrid.Agent(:,:,2);          
voxelgrid.Agent(end-1,:,:) = voxelgrid.Agent(end,:,:)|voxelgrid.Agent(end-1,:,:);          
voxelgrid.Agent(:,end-1,:) = voxelgrid.Agent(:,end,:)|voxelgrid.Agent(:,end-1,:);          
voxelgrid.Agent(:,:,end-1) = voxelgrid.Agent(:,:,end)|voxelgrid.Agent(:,:,end-1);

figure
hold on
        xlim([0 500])%Set lims so matlab does not graph over Mv
        ylim([0 500])
        zlim([0 500])
p=patch(isosurface(voxelgrid.Agent==1,0));
set(p,'facecolor','red' ,'edgecolor', 'none');
daspect([1 1 1])
%rotate3d on;
view([60 30]);
camlight
lighting gouraud



CoolPlotRadius = 1; % allows smaller than 1 so can see inside...
Coolcellradius = 5;
[spherex,spherey,spherez] = sphere;% (8);
spherex = spherex.*CoolPlotRadius.*Coolcellradius;
spherey = spherey.*CoolPlotRadius.*Coolcellradius;
spherez = spherez.*CoolPlotRadius.*Coolcellradius;
hold all;

disp(strcat("plotting ",num2str(length(XYZ))," cells"));

XYZ=10*XYZ;
scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'filled','MarkerEdgeColor','k')

disp(strcat("plotting ",num2str(length(XYZkt))," kt cells"));
XYZkt2 = XYZkt(:,:);
XYZkt2 = XYZkt2*10;
scatter3(XYZkt2(:,1),XYZkt2(:,2),XYZkt2(:,3),'filled','MarkerEdgeColor','k')

%% Adjust axes settings
%brighten(0.2);
light
%lighting phong;
box on;		
grid off;

axis equal; 
%set(gca,'position',[0 0 1 1]);

rotate3d on;		
hold off;		
axis equal; 
axis([0 500, 0 500, 0 500]); %set 2nd 4th and 6th vector elements to manually fix plot size
view(150,36);

end


% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2
value1=get(hObject, 'Value');

if(value1==1)
    disp('2!');

else
    disp('nothing again');
end