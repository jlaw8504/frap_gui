function varargout = FRAP_gui(varargin)
% FRAP_GUI MATLAB code for FRAP_gui.fig
%      FRAP_GUI, by itself, creates a new FRAP_GUI or raises the existing
%      singleton*.
%
%      H = FRAP_GUI returns the handle to a new FRAP_GUI or the handle to
%      the existing singleton*.
%
%      FRAP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRAP_GUI.M with the given input arguments.
%
%      FRAP_GUI('Property','Value',...) creates a new FRAP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FRAP_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FRAP_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FRAP_gui

% Last Modified by GUIDE v2.5 18-Jan-2017 11:54:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FRAP_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @FRAP_gui_OutputFcn, ...
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


% --- Executes just before FRAP_gui is made visible.
function FRAP_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FRAP_gui (see VARARGIN)

% Choose default command line output for FRAP_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FRAP_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FRAP_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pre_button.
function pre_button_Callback(hObject, eventdata, handles)
% hObject    handle to pre_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Open Pre-laser Image Stack
pre_file = uigetfile('*.*');
pre_data = bfopen(pre_file);
%% Isolate GFP (Color 1)
%parse the plane, Z, and color information from data
im_info = pre_data{1,1}(:,2);
%find GFP image index
gfp_idx = cellfun(@(x) ~isempty(strfind(x,'C=1/3')),im_info);
pre_gfp = pre_data{1,1}(gfp_idx,1);
handles.pre_gfp = pre_gfp;

%% setup for pre_slider
%get plane number
pre_plane_num = size(pre_gfp,1);
handles.pre_plane_num = pre_plane_num;
set(handles.pre_slider, 'SliderStep', [1/(pre_plane_num-1) 1/(pre_plane_num-1)],...
    'Min', 1, 'Max', pre_plane_num, 'Value', ceil(pre_plane_num/2));
%% Show first image of the stack
%parse slider Value
pre_slider_value = get(handles.pre_slider,'Value');
%plot image in axes1
subplot(handles.axes1);
imshow(pre_gfp{pre_slider_value},[]);
%upate handles structure
guidata(hObject,handles);


% --- Executes on button press in laser_button.
function laser_button_Callback(hObject, eventdata, handles)
% hObject    handle to laser_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Open Laser Image 
laser_file = uigetfile('*.*');
laser_data = bfopen(laser_file);
%% Isolate GFP (Color 2)
%parse the plane, Z, and color information from data
im_info = laser_data{1,1}(:,2);
%find GFP image index
gfp_idx = cellfun(@(x) ~isempty(strfind(x,'C=2/2')),im_info);
laser_gfp = laser_data{1,1}{gfp_idx,1};
handles.laser_gfp = laser_gfp;

%% Show image
%plot image in axes2
subplot(handles.axes2);
imshow(laser_gfp,[]);
%upate handles structure
guidata(hObject,handles);

% --- Executes on button press in post_button.
function post_button_Callback(hObject, eventdata, handles)
% hObject    handle to post_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Open Post-laser Image Stack
post_file = uigetfile('*.*');
post_data = bfopen(post_file);
post_gfp = post_data{1,1}(:,1);
handles.post_gfp = post_gfp;

%% setup for post_slider
%get plane number
post_plane_num = size(post_gfp,1);
handles.post_plane_num = post_plane_num;
set(handles.post_slider, 'SliderStep', [1/(post_plane_num-1) 1/(post_plane_num-1)],...
    'Min', 1, 'Max', post_plane_num, 'Value', 1);
%% Show first image of the stack
%parse slider Value
post_slider_value = get(handles.post_slider,'Value');
%plot image in axes3
subplot(handles.axes3);
imshow(post_gfp{post_slider_value},[]);
%upate handles structure
guidata(hObject,handles);

% --- Executes on slider movement.
function pre_slider_Callback(hObject, eventdata, handles)
% hObject    handle to pre_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%% Update axes1 with slider values
pre_slider_value = get(handles.pre_slider,'Value');
%imshow to axes1
subplot(handles.axes1);
imshow(handles.pre_gfp{pre_slider_value},[]);
%upate handles structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pre_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pre_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function post_slider_Callback(hObject, eventdata, handles)
% hObject    handle to post_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%% Update axes1 with slider values
post_slider_value = get(handles.post_slider,'Value');
%imshow to axes3
subplot(handles.axes3);
imshow(handles.post_gfp{post_slider_value},[]);
%upate handles structure
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function post_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to post_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in bin_button.
function bin_button_Callback(hObject, eventdata, handles)
% hObject    handle to bin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Threshold to create binary mask of laser image
%Threshold
thresh = multithresh(handles.laser_gfp);
%make binary mask
im_pre_bin = handles.laser_gfp > thresh;
%set all zero values in im_pre_bin to min values of im_pre_bin index
im_parsed = handles.laser_gfp(im_pre_bin);
im_min = min(im_parsed(:));
handles.laser_gfp(~im_pre_bin) = im_min;
%repeat thresholding to generate im_bin
thresh2 = multithresh(handles.laser_gfp);
im_bin = handles.laser_gfp > thresh2;
%plot binary mask in axes2
subplot(handles.axes2);
imshow(im_bin);
%add im_bin to handles
handles.im_bin = im_bin;
%upate handles structure
guidata(hObject,handles);


% --- Executes on button press in sig_bg_button.
function sig_bg_button_Callback(hObject, eventdata, handles)
% hObject    handle to sig_bg_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number for pre_slider
slider_value = get(handles.pre_slider,'Value');
%change the plot to axes1
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle
pos_sig_bg = ceil(wait(h));
handles.pos_sig_bg = pos_sig_bg;
%sample the GFP image
im_sig_bg = handles.pre_gfp{slider_value}(round(pos_sig_bg(2)):round(pos_sig_bg(2))+round((pos_sig_bg(4))),round(pos_sig_bg(1)):...
    round(pos_sig_bg(1))+ round(pos_sig_bg(3)));
%get average of that image
sig_bg = round(mean(im_sig_bg(:)));
handles.sig_bg = sig_bg;
set(handles.sig_bg_text, 'String', num2str(sig_bg));
guidata(hObject,handles);

% --- Executes on button press in ref_button.
function ref_button_Callback(hObject, eventdata, handles)
% hObject    handle to ref_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number for pre_slider
slider_value = get(handles.pre_slider,'Value');
%change the plot to axes1
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle
pos_ref = ceil(wait(h));
handles.pos_ref = pos_ref;
%sample the GFP image
im_ref = handles.pre_gfp{slider_value}(round(pos_ref(2)):round(pos_ref(2))+round((pos_ref(4))),round(pos_ref(1)):...
    round(pos_ref(1))+ round(pos_ref(3)));
%get average of that image
ref = round(mean(im_ref(:)));
handles.ref = ref;
set(handles.ref_text, 'String', num2str(ref));
guidata(hObject,handles);

% --- Executes on button press in ref_bg_button.
function ref_bg_button_Callback(hObject, eventdata, handles)
% hObject    handle to ref_bg_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number for pre_slider
slider_value = get(handles.pre_slider,'Value');
%change the plot to axes1
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle
pos_ref_bg = ceil(wait(h));
handles.pos_ref_bg = pos_ref_bg;
%sample the GFP image
im_ref_bg = handles.pre_gfp{slider_value}(round(pos_ref_bg(2)):round(pos_ref_bg(2))+round((pos_ref_bg(4))),round(pos_ref_bg(1)):...
    round(pos_ref_bg(1))+ round(pos_ref_bg(3)));
%get average of that image
ref_bg = round(mean(im_ref_bg(:)));
handles.ref_bg = ref_bg;
set(handles.ref_bg_text, 'String', num2str(ref_bg));
guidata(hObject,handles);


% --- Executes on button press in sig_button.
function sig_button_Callback(hObject, eventdata, handles)
% hObject    handle to sig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number for pre_slider
slider_value = get(handles.pre_slider,'Value');
%change the plot to axes1
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle
pos_sig = ceil(wait(h));
handles.pos_sig = pos_sig;
%sample the GFP image
im_sig = handles.pre_gfp{slider_value}(round(pos_sig(2)):round(pos_sig(2))+round((pos_sig(4))),round(pos_sig(1)):...
    round(pos_sig(1))+ round(pos_sig(3)));
%get average of that image
sig = round(mean(im_sig(:)));
handles.sig = sig;
set(handles.sig_text, 'String', num2str(sig));
guidata(hObject,handles);


% --- Executes on button press in sub_bin_button.
function sub_bin_button_Callback(hObject, eventdata, handles)
% hObject    handle to sub_bin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add binary mask to pre image
%load binary mask
im_bin = handles.im_bin;
%load pre_slider value
pre_slider_value = get(handles.pre_slider,'Value');
%load correct gfp image plane
im = handles.pre_gfp{pre_slider_value};
%find max of image
im_max = max(im(:));
%use logical indexing to set im_bin pixels to max intensity
im(im_bin) = im_max;
%update axes1 to show new image
subplot(handles.axes1);
imshow(im,[]);


% --- Executes on button press in frap_button.
function frap_button_Callback(hObject, eventdata, handles)
% hObject    handle to frap_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Calculate FRAP curves
%load binary masks
im_bin = handles.im_bin;
im_bin_resize = handles.im_bin_resize;
%create mask of unbleached portion or nucleus
im_unbleach_bin = im_bin + im_bin_resize;
im_unbleach_bin = im_unbleach_bin == 1;

%load all positional information
pos_sig_bg = handles.pos_sig_bg;
pos_ref = handles.pos_ref;
pos_ref_bg = handles.pos_ref_bg;
%% Calculate mean of all positions in pre_laser and post_laser images
total_planes = handles.post_plane_num + 1;
%pre-allocate vector
mean_sig = zeros([total_planes 1]);
mean_sig_bg = zeros([total_planes 1]);
mean_unbleach = zeros([total_planes 1]);
mean_ref = zeros([total_planes 1]);
mean_ref_bg = zeros([total_planes 1]);
for n = 1:(total_planes)
    if n == 1
%load in correct image
im = handles.pre_gfp{get(handles.pre_slider,'Value')};
    else
        im = handles.post_gfp{n-1};
    end
%load in positions and parse pre_im
sig_bg = im(round(pos_sig_bg(2)):round(pos_sig_bg(2))+round((pos_sig_bg(4))),round(pos_sig_bg(1)):...
    round(pos_sig_bg(1))+ round(pos_sig_bg(3)));
ref = im(round(pos_ref(2)):round(pos_ref(2))+round((pos_ref(4))),round(pos_ref(1)):...
    round(pos_ref(1))+ round(pos_ref(3)));
ref_bg = im(round(pos_ref_bg(2)):round(pos_ref_bg(2))+round((pos_ref_bg(4))),round(pos_ref_bg(1)):...
    round(pos_ref_bg(1))+ round(pos_ref_bg(3)));
%calc mean intensity value of each parsed position
mean_sig(n,1) = mean(im(im_bin));
mean_sig_bg(n,1) = mean(sig_bg(:));
mean_unbleach(n,1) = mean(im(im_unbleach_bin));
mean_ref(n,1) = mean(ref(:));
mean_ref_bg(n,1) = mean(ref_bg(:));
end

%% Calculate BG corr values
mean_sig_bg_corr = mean_sig - mean_sig_bg;
mean_unbleach_bg_corr = mean_unbleach - mean_sig_bg;
mean_ref_bg_corr = mean_ref - mean_ref_bg;
%% Photobleach correction
%since first set is image stack only correct bleaching in post timelapse
timepoints = (1:size(mean_ref_bg_corr,1)-1)';
linear_coeff = polyfit(timepoints,mean_ref_bg_corr(2:end),1);
ref_slope = linear_coeff(1);
%generate correction vector for photobleaching
%skip over first entry (pre-bleach stack will not be corrected)
corr_vect = -ref_slope .* timepoints;
corr_vect = [0;corr_vect];
sig_bg_photo_corr = mean_sig_bg_corr + corr_vect;
unbleach_bg_photo_corr = mean_unbleach_bg_corr + corr_vect;
ref_bg_photo_corr = mean_ref_bg_corr + corr_vect;
%% Create normalized signal and reference data
sig_norm = sig_bg_photo_corr./sig_bg_photo_corr(1);
unbleach_norm = unbleach_bg_photo_corr./unbleach_bg_photo_corr(1);
ref_norm =ref_bg_photo_corr./ref_bg_photo_corr(1);
%plot out the normalized data
subplot(handles.axes4);
plot(sig_norm);
hold on;
plot(unbleach_norm);
plot(ref_norm);
legend('Bleached', 'Unbleached', 'Reference');
hold off;
%% Save data to MAT file
%load in relevant data
bleach_fraction = handles.bleach_fraction;
nuc_area = handles.nuc_area;
laser_area = handles.laser_area;
save('FRAP.mat', 'pos*', 'sig*', 'unbleach*',...
    '*ref*', 'bleach_fraction', '*_area');


% --- Executes on button press in erode_button.
function erode_button_Callback(hObject, eventdata, handles)
% hObject    handle to erode_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Erode im_bin
%use disk structural element
se_disk = strel('disk',1');
handles.im_nuc_bin = imerode(handles.im_nuc_bin,se_disk);
subplot(handles.axes5);
imshow(handles.im_nuc_bin,[]);
%Update handles
guidata(hObject,handles);


% --- Executes on button press in nuclear_button.
function nuclear_button_Callback(hObject, eventdata, handles)
% hObject    handle to nuclear_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the current plane number for pre_slider
slider_value = get(handles.pre_slider,'Value');
%change the plot to axes1
subplot(handles.axes1);
%draw a rectangle
h = imrect;
%get position coordinates of rectangle and add to handles
pos_nuc = ceil(wait(h));
handles.pos_nuc = pos_nuc;
%sample the GFP image
im_nuc = handles.pre_gfp{slider_value}(round(pos_nuc(2)):round(pos_nuc(2))+round((pos_nuc(4))),round(pos_nuc(1)):...
    round(pos_nuc(1))+ round(pos_nuc(3)));
handles.im_nuc = im_nuc;
mean_im_nuc = round(mean(im_nuc(:)));
handles.mean_im_nuc = mean_im_nuc;
set(handles.nuclear_text, 'String', num2str(mean_im_nuc));
guidata(hObject,handles);


% --- Executes on button press in nuc_bin_button.
function nuc_bin_button_Callback(hObject, eventdata, handles)
% hObject    handle to nuc_bin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Use BG sub and thresholding to generate a nuclear binary for area
%BG sub
im_nuc_sub = handles.im_nuc - handles.sig_bg;
%threshold using multithresh
thresh = multithresh(im_nuc_sub);
%create nuclear binary image
im_nuc_bin = im_nuc_sub > thresh;
%display nuclear binary on axes5 for user inspection and alteration
subplot(handles.axes5);
imshow(im_nuc_bin);
handles.im_nuc_bin = im_nuc_bin;
%Update handles structure
guidata(hObject,handles);

% --- Executes on button press in dilate_button.
function dilate_button_Callback(hObject, eventdata, handles)
% hObject    handle to dilate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Dilate im_bin
%use disk structural element
se_disk = strel('disk',1');
handles.im_nuc_bin = imdilate(handles.im_nuc_bin,se_disk);
subplot(handles.axes5);
imshow(handles.im_nuc_bin,[]);
%Update handles
guidata(hObject,handles);


% --- Executes on button press in overlay_nuc_button.
function overlay_nuc_button_Callback(hObject, eventdata, handles)
% hObject    handle to overlay_nuc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Add binary mask to pre image
%load in im_nuc_bin and resize to original image
im_nuc_bin = handles.im_nuc_bin;
%get pos_nuc
pos_nuc = handles.pos_nuc;
%load pre_slider value
pre_slider_value = get(handles.pre_slider,'Value');
%load correct gfp image plane
im = handles.pre_gfp{pre_slider_value};
%% Pad the binary image back to the orginal size
[im_y, im_x, ~] = size(im);
x_pad_1 = zeros(pos_nuc(4) + 1, pos_nuc(1) - 1);
im_bin_resize = horzcat(x_pad_1,im_nuc_bin);
[~,im_bin_resize_width] = size(im_bin_resize);
x_pad_2 = zeros(pos_nuc(4) +1, im_x - im_bin_resize_width);
im_bin_resize = horzcat(im_bin_resize, x_pad_2);

y_pad_1 = zeros(pos_nuc(2) - 1, im_x);
im_bin_resize = vertcat(y_pad_1,im_bin_resize);
[im_bin_resize_height, ~] = size(im_bin_resize);
y_pad_2 = zeros(im_y - im_bin_resize_height, im_x);
im_bin_resize = vertcat(im_bin_resize, y_pad_2);

%find max of image
im_max = max(im(:));
%use logical indexing to set im_bin pixels to max intensity
im(logical(im_bin_resize)) = im_max;
%update axes1 to show new image
subplot(handles.axes1);
imshow(im,[]);
%Add im_bin_resize to handles and update handles
handles.im_bin_resize = im_bin_resize;
guidata(hObject,handles);

% --- Executes on button press in area_button.
function area_button_Callback(hObject, eventdata, handles)
% hObject    handle to area_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Calculate percent of nucleus bleached
%calculate areas of binary images
nuc_stats = regionprops(handles.im_nuc_bin,'Area');
laser_stats = regionprops(handles.im_bin,'Area');
bleach_fraction = laser_stats.Area/nuc_stats.Area;
set(handles.area_text, 'String', num2str(round(bleach_fraction,2)));
%Update handles
handles.nuc_area = nuc_stats.Area;
handles.laser_area = laser_stats.Area;
handles.bleach_fraction = bleach_fraction;
guidata(hObject,handles);
