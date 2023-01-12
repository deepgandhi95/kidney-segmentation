function varargout = MRE_KIDNEY_DEEP(varargin)
% MRE_KIDNEY_DEEP MATLAB code for MRE_KIDNEY_DEEP.fig
%      MRE_KIDNEY_DEEP, by itself, creates a new MRE_KIDNEY_DEEP or raises the existing
%      singleton*.
%
%      H = MRE_KIDNEY_DEEP returns the handle to a new MRE_KIDNEY_DEEP or the handle to
%      the existing singleton*.
%
%      MRE_KIDNEY_DEEP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRE_KIDNEY_DEEP.M with the given input arguments.
%
%      MRE_KIDNEY_DEEP('Property','Value',...) creates a new MRE_KIDNEY_DEEP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MRE_KIDNEY_DEEP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MRE_KIDNEY_DEEP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MRE_KIDNEY_DEEP

% Last Modified by GUIDE v2.5 09-Feb-2017 16:20:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRE_KIDNEY_DEEP_OpeningFcn, ...
                   'gui_OutputFcn',  @MRE_KIDNEY_DEEP_OutputFcn, ...
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


% --- Executes just before MRE_KIDNEY_DEEP is made visible.
function MRE_KIDNEY_DEEP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRE_KIDNEY_DEEP (see VARARGIN)

% Choose default command line output for MRE_KIDNEY_DEEP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MRE_KIDNEY_DEEP wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Begin the pointer with an arrow
set(gcf, 'pointer', 'arrow')

% --- Outputs from this function are returned to the command line.
function varargout = MRE_KIDNEY_DEEP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open_localizer.
function open_localizer_Callback(hObject, eventdata, handles)
% hObject    handle to open_localizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_array;
global index; index = 1; 
global filenames; global pathstr;  
global data; 
global size_offset; global size_slice; global st; 
global magX; global magY;  
global phaseX; global phaseY; global phaseZ; 
global masksX; global masksY; global masksZ;  

set(gcf, 'pointer', 'watch')
drawnow;

[filenames, pathstr] = uigetfile('*.*', 'Select MRE data or loc file','multiselect', 'on');
% end
tmp = dicomread(fullfile(pathstr, filenames{1,1})); 

data = zeros(size(tmp,1),size(tmp,2),size(filenames,2)); 
for i = 1:size(filenames,2)
    
data(:,:,i) = dicomread(fullfile(pathstr, filenames{1,i}));
end
 
imagesc(data(:,:,index)); colormap(gray);axis off; hold on;
title(sprintf('Magnitude Z index = %.3g%', index), 'fontsize', 18);
set(gcf, 'pointer', 'arrow')

masksZ = zeros(size(tmp,1), size(tmp,2), 2,size(filenames,2)); 


function keypress_callback(hObject, eventdata, handles)

global index; global data; 
global mask1; global mask2; 
global filenames; 
tmp = filenames; 
global masksZ
%%%%For normal order of files%%%%

    
switch eventdata.Character
    case 28 % Left
        if index == 1
            index = size(tmp,3); 
        else 
            index = index - 1; 
        end 
    case 29 % Right
        if index == size(tmp,2)
            index = 1; 
        else 
            index = index + 1; 
        end 
end 


imagesc(data(:,:,index)); colormap(gray);axis off; hold on;
title(sprintf('Magnitude Z index = %.3g%', index), 'fontsize', 18); 

mask1 = masksZ(:,:,1,index); 
mask2 = masksZ(:,:,2,index); 


    if max(max(mask1))~=0
    contour(mask1,'r','linewidth',0.5);hold on; 
    contour(mask2,'r','linewidth',0.5);hold off; 
    else 
        hold off; 
    end




% --- Executes during object creation, after setting all properties.
function biggest_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to biggest_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

% --- Executes on button press in open_vars.
function open_vars_Callback(hObject, eventdata, handles)
% hObject    handle to open_vars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global file_array;
global index; index = 1; 
global filenames; global pathstr;  
global data1; 
global size_offset; global size_slice; global st; 
global magX; global magY; global magZ; 
global phaseX; global phaseY; global phaseZ; 
global masksX; global masksY; global masksZ;  

%Change the pointer to watch when the vars file is loading
set(gcf, 'pointer', 'watch')
drawnow;



% % %Select vars data and load it. Load last directory, if present. Save to temp
% % file.
% if exist('\\Svr-dfskrc-vp01\dfs\Share3\Shared.Services\RadiologyLLC\Kolipaka-Lab\temp\prateek\temp\lstdr5.mat','file')
%     load('\\Svr-dfskrc-vp01\dfs\Share3\Shared.Services\RadiologyLLC\Kolipaka-Lab\temp\prateek\temp\lstdr5.mat');
%     [filenames,pathstr] = uigetfile(fullfile(pname,'*.*'),'Select the MRE data', 'MultiSelect', 'on');
% else
[filenames, pathstr] = uigetfile('*.*', 'Select MRE data or vars file','multiselect', 'on');
% end

%Select and load vars file
% [filenames, pathstr] = uigetfile('*.*', 'Select MRE data or vars file','multiselect', 'on');

data1 = load(fullfile(pathstr, filenames)); 
file_array = fieldnames(data1);
file_array = sort(file_array); 

% --- Executes on button press in load_mask.
function load_mask_Callback(hObject, eventdata, handles)
% hObject    handle to load_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global masksX; global masksY; global masksZ;  

[file_name,p,~] = uigetfile('*.*');
load(fullfile(p,file_name)); 



% --- Executes on button press in slice.
function slice_Callback(hObject, eventdata, handles)
% hObject    handle to slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slice


% --- Executes on button press in Slice_and_Offsets.
function Slice_and_Offsets_Callback(hObject, eventdata, handles)
% hObject    handle to Slice_and_Offsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Slice_and_Offsets


% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global index;  global data; 
global mask1; global mask2;  
global beginslice; global masksZ; 


beginslice = 1;%str2double(get(handles.biggest_slice,'String')); 

imagesc(data(:,:,index)); colormap(gray);axis off; hold on;
title(sprintf('Magnitude Z index = %.3g%', index), 'fontsize', 18);
mask1 = masksZ(:,:,1,index); 
mask2 = masksZ(:,:,2,index); 


if max(max(mask1))~=0
contour(mask1,'r','linewidth',0.5);hold on; 
contour(mask2,'r','linewidth',0.5);hold off; 
else 
hold off; 
end


set(gcf, 'KeyPressFcn', @keypress_callback, 'Interruptible', 'on', 'BusyAction', 'queue'); 


% --- Executes on button press in select_roi.
function select_roi_Callback(hObject, eventdata, handles)
% hObject    handle to select_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global magX; global magY; global magZ; global index; 
global data;
global roi1; global roi2; 


tmpim = data(:,:,index); 
 
figure, imshow(tmpim, []); 
title('Select first ROI: Left Kidney'); 

roi1 = imfreehand; 
roi1= createMask(roi1); 
close(gcf); 
 

figure, imshow(tmpim, []); 
title('Select second ROI: Right Kidney'); 

roi2 = imfreehand; 
roi2 = createMask(roi2); 
close(gcf); 

% --- Executes on button press in erosion.
function erosion_Callback(hObject, eventdata, handles)
% hObject    handle to erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
global mask1; global mask2; 
global masksX; global masksY; global masksZ; 
global size_offset; 
global data;
% Construct a questdlg 
choice = questdlg('This will erode the boundary pixels on selected direction. Do you still want to do it?', ...
'STOP', ...
'Yes','No','No');
% Handle response      
 
   mask1 = data(:,:,1,index); 
   mask2 = data(:,:,2,index); 
      
    E = bwmorph(mask1, 'endpoints'); 
    E = find(E==1);
    mask1(E) = 0; 

    E = bwmorph(mask2, 'endpoints'); 
    E = find(E==1);
    mask2(E) = 0; 
    
    %%Pick biggest area (This would eliminate the small un-connected
    %%regions that are created as a result of erosion
    s1 = mask1; 
    cc = bwconncomp(s1);
    stats = regionprops(cc, 'basic');
    A = [stats.Area];
    [~, biggest] = max(A);
    s1(labelmatrix(cc)~=biggest) = 0; 
    s1 = imfill(s1,'holes');
    mask1 = s1; 

    %%Pick biggest area
    s1 = mask2; 
    cc = bwconncomp(s1);
    stats = regionprops(cc, 'basic');
    A = [stats.Area];
    [~, biggest] = max(A);
    s1(labelmatrix(cc)~=biggest) = 0; 
    s1 = imfill(s1,'holes');
    mask2 = s1; 
    
 data(:,:,1,i)=mask1; 
 data(:,:,2,i)=mask2; 
    
    %Do nothing
update_mask(handles); 


% --- Executes on button press in save_masked_vars_and_masks.
function save_masked_vars_and_masks_Callback(hObject, eventdata, handles)
% hObject    handle to save_masked_vars_and_masks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data1; global masksZ; global file_array; global pathstr;
mask_l = squeeze(masksZ(:,:,1,:));
mask_r = squeeze(masksZ(:,:,2,:));
mask_l = cat(4,mask_l,mask_l,mask_l,mask_l);
mask_l = mask_l | mask_l;
mask_r = cat(4,mask_r,mask_r,mask_r,mask_r);
mask_r = mask_r | mask_r;
mask_lr = mask_l | mask_r;


xx = data1.(file_array{1});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','M'));
eval([genvarname '= tmp;']);
xx = data1.(file_array{2});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','P'));
eval([genvarname '= tmp;']);
xx = data1.(file_array{3});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','M'));
eval([genvarname '= tmp;']);
xx = data1.(file_array{4});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','P'));
eval([genvarname '= tmp;']);
xx = data1.(file_array{5});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','M'));
eval([genvarname '= tmp;']);
xx = data1.(file_array{6});
tmp = mask_l(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','P'));
eval([genvarname '= tmp;']);

filename=char(strcat(pathstr,'\vars_masked_l.mat'));
save(filename, 'DIR_*');%save vars_masked file that contains all files begin with DIR_....
filename_masks=char(strcat(pathstr,'\masks_l.mat'));%save masks file
save(filename_masks, 'masks*');

clearvars('DIR_X_M','DIR_X_P','DIR_Y_M','DIR_Y_P','DIR_Z_M','DIR_Z_P');

xx = data1.(file_array{1});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','M'));
eval([genvarname '= tmp1;']);
xx = data1.(file_array{2});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','P'));
eval([genvarname '= tmp1;']);
xx = data1.(file_array{3});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','M'));
eval([genvarname '= tmp1;']);
xx = data1.(file_array{4});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','P'));
eval([genvarname '= tmp1;']);
xx = data1.(file_array{5});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','M'));
eval([genvarname '= tmp1;']);
xx = data1.(file_array{6});
tmp1 = mask_r(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','P'));
eval([genvarname '= tmp1;']);

filename=char(strcat(pathstr,'\vars_masked_r.mat'));
save(filename, 'DIR_*');%save vars_masked file that contains all files begin with DIR_....
filename_masks=char(strcat(pathstr,'\masks_r.mat'));%save masks file
save(filename_masks, 'masks*');

clearvars('DIR_X_M','DIR_X_P','DIR_Y_M','DIR_Y_P','DIR_Z_M','DIR_Z_P');

xx = data1.(file_array{1});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','M'));
eval([genvarname '= tmp2;']);
xx = data1.(file_array{2});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','X','_','P'));
eval([genvarname '= tmp2;']);
xx = data1.(file_array{3});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','M'));
eval([genvarname '= tmp2;']);
xx = data1.(file_array{4});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Y','_','P'));
eval([genvarname '= tmp2;']);
xx = data1.(file_array{5});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','M'));
eval([genvarname '= tmp2;']);
xx = data1.(file_array{6});
tmp2 = mask_lr(:,:,:,:).*xx;
genvarname=char(strcat('DIR_','Z','_','P'));
eval([genvarname '= tmp2;']);

filename=char(strcat(pathstr,'\vars_masked_lr.mat'));
save(filename, 'DIR_*');%save vars_masked file that contains all files begin with DIR_....
filename_masks=char(strcat(pathstr,'\masks_lr.mat'));%save masks file
save(filename_masks, 'masks*');

clearvars('DIR_X_M','DIR_X_P','DIR_Y_M','DIR_Y_P','DIR_Z_M','DIR_Z_P');

pname = pathstr; 

% --- Executes on button press in add_1.
function add_1_Callback(hObject, eventdata, handles)
% hObject    handle to add_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global index; global mask1;
global magX; global magY; global magZ; 
global masksX; global masksY; global masksZ; 
global data;

tmpim = data(:,:,index); 
mask = masksZ(:,:,1,index);


cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask1,'r','linewidth',0.5);hold off; 
axes(handles.axes1);  
h = imfreehand;
bw = createMask(h); 
mask = mask | bw; 
delete(h); 
mask1 = mask; 

masksZ(:,:,1,index) = mask1; 
  
update_mask(handles); 

% --- Executes on button press in remove_1.
function remove_1_Callback(hObject, eventdata, handles)
% hObject    handle to remove_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; global mask1;
global magX; global magY; global magZ; 
global data;
global masksX; global masksY; global masksZ; 



tmpim = data(:,:,index); 
mask = masksZ(:,:,1,index);


cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask1,'r','linewidth',0.5);hold off; 
axes(handles.axes1);  
h = imfreehand;
bw = createMask(h); 
mask(bw) = 0; 
delete(h); 
mask1 = mask; 

masksZ(:,:,1,index) = mask1; 

update_mask(handles); 


% --- Executes on button press in biggest_1.
function biggest_1_Callback(hObject, eventdata, handles)
% hObject    handle to biggest_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; global mask1;
global magX; global magY; global magZ; 
global data;
global masksX; global masksY; global masksZ; 


tmpim = data(:,:,index); 

cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask1,'r','linewidth',0.5);hold off; 
%%Pick biggest area
s1 = mask1; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask1 = s1; 

masksZ(:,:,1,index) = mask1; 
 
update_mask(handles); 

% --- Executes on button press in Threshold1.
function Threshold1_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; 
global data; 
global roi1;
global mask1; 
global level; 
global masksZ; 


tmpim = data(:,:,index); 


a = tmpim; 
figure, imshow(a,[]); 
level = thresh_tool(a, []);  

max_val = max(max(a)); 

img = a; 

tmp = zeros(size(img));
e = find( (img>level) );
tmp(e) = max_val; 

% roi1 = imfreehand; 
% roi1= createMask(roi1); 
% close(gcf);
mask1 = tmp.*roi1; 

%%Pick biggest area
s1 = mask1; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask1 = s1;

%Save masks and update

masksZ(:,:,1,index) = mask1; 

update_mask(handles); 

% --- Executes on button press in add_2.
function add_2_Callback(hObject, eventdata, handles)
% hObject    handle to add_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; global mask2;
global magX; global magY; global magZ; 
global data;
global masksX; global masksY; global masksZ; 


tmpim = data(:,:,index); 
mask = masksZ(:,:,2,index);

cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask2,'r','linewidth',0.5);hold off; 
axes(handles.axes1);  
h = imfreehand;
bw = createMask(h); 
mask = mask | bw; 
delete(h); 
mask2 = mask; 

masksZ(:,:,2,index) = mask2; 
 
update_mask(handles); 


% --- Executes on button press in remove_2.
function remove_2_Callback(hObject, eventdata, handles)
% hObject    handle to remove_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; global mask2;
global magX; global magY; global magZ; 
global data;
global masksX; global masksY; global masksZ; 



tmpim = data(:,:,index); 
mask = masksZ(:,:,2,index);


cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask2,'r','linewidth',0.5);hold off; 
axes(handles.axes1);  
h = imfreehand;
bw = createMask(h); 
mask(bw) = 0; 
delete(h); 
mask2 = mask; 

masksZ(:,:,2,index) = mask2; 

update_mask(handles); 

% --- Executes on button press in biggest_2.
function biggest_2_Callback(hObject, eventdata, handles)
% hObject    handle to biggest_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global index; global mask2;
global magX; global magY; global magZ; 
global data;
global masksX; global masksY; global masksZ; 



tmpim = data(:,:,index); 

cIM = tmpim; 
cIM = uint16(cIM);
axes(handles.axes1);  
imagesc(cIM); colormap(gray);axis off; hold on; 
contour(mask2,'r','linewidth',0.5);hold off; 
%%Pick biggest area
s1 = mask2; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask2 = s1; 

masksZ(:,:,2,index) = mask2; 

update_mask(handles);

% --- Executes on button press in Threshold_2.
function Threshold_2_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global magX; global magY; global magZ; global index; 
global data
global roi2; 
global mask2; 
global masksX; global masksY; global masksZ; 

tmpim = data(:,:,index); 

a = tmpim; 
figure, imshow(a,[]); 
level = thresh_tool(a, []);  

max_val = max(max(a)); 

img = a; 

tmp = zeros(size(img));
e = find( (img>level) );
tmp(e) = max_val; 

% roi2 = imfreehand; 
% roi2 = createMask(roi2); 
% close(gcf);
mask2 = tmp.*roi2; 

%%Pick biggest area
s1 = mask2; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask2 = s1;

%Save masks and update

masksZ(:,:,2,index) = mask2; 

update_mask(handles);

% --- Executes on button press in magX.
function magX_Callback(hObject, eventdata, handles)
% hObject    handle to magX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of magX


% --- Executes on button press in Mag_Y.
function Mag_Y_Callback(hObject, eventdata, handles)
% hObject    handle to Mag_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mag_Y


% --- Executes on button press in mag_z.
function mag_z_Callback(hObject, eventdata, handles)
% hObject    handle to mag_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mag_z

function update_mask(handles)
%Show masked boundaries on top of image

global magX; global magY; global magZ; global index; 
global mask1; global mask2;
 global masksZ; global data; 


tmpim = data(:,:,index); 
mask1 = masksZ(:,:,1,index); 
mask2 = masksZ(:,:,2,index); 

axes(handles.axes1);

imagesc(tmpim); colormap(gray);axis off; hold on; 
if max(max(mask1))~=0
    contour(mask1,'r','linewidth',0.5);hold on; 
end
if max(max(mask2))~=0
    contour(mask2,'r','linewidth',0.5);hold off; 
else
    hold off;
end                                                                                                                                         


% --- Executes on button press in Threshold_all.
function Threshold_all_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global magX; global magY; global magZ;
global filenames;
global roi1; global roi2; 
global mask1; global mask2; 
global level; 
global masksZ; 
global data; global index;

for i = 1:size(filenames,2)

tmpim = data(:,:,i); 

a = tmpim; 
max_val = max(max(a)); 

img = a; 

tmp = zeros(size(img));
e = find( (img>level) );
tmp(e) = max_val; 

mask1 = tmp.*roi1; 

%%Pick biggest area
s1 = mask1; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask1 = s1;

mask2 = tmp.*roi2; 

%%Pick biggest area
s1 = mask2; 
cc = bwconncomp(s1);
stats = regionprops(cc, 'basic');
A = [stats.Area];
[~, biggest] = max(A);
s1(labelmatrix(cc)~=biggest) = 0; 
s1 = imfill(s1,'holes');
mask2 = s1;

%save maskes to only only selected direction

masksZ(:,:,1,i)=mask1; 
masksZ(:,:,2,i)=mask2; 

end

update_mask(handles); 