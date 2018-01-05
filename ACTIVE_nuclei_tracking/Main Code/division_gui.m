function varargout = division_gui(varargin)
% DIVISION_GUI MATLAB code for division_gui.fig
%
% To run: 
% [div_output, pause_status] = division_gui(filename, event_array, xyzs_id, xyzs_id_columns);
%
% To run previously saved file:
% [div_output, pause_status] = division_gui(filename, event_array, xyzs_id, xyzs_id_columns, currentEventIterator, currentEventNum, div_output);
%
% Megan Brasch, Richard Baker  12/9/13
%
%      DIVISION_GUI, by itself, creates a new DIVISION_GUI or raises the existing
%      singleton*.
%
%      H = DIVISION_GUI returns the handle to a new DIVISION_GUI or the handle to
%      the existing singleton*.
%
%      DIVISION_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIVISION_GUI.M with the given input arguments.
%
%      DIVISION_GUI('Property','Value',...) creates a new DIVISION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before division_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to division_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help division_gui

% Last Modified by GUIDE v2.5 09-Jan-2014 10:34:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @division_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @division_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
warning('off','MATLAB:str2func:invalidFunctionName')

if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before division_gui is made visible.
function division_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to division_gui (see VARARGIN)

% Allow relevant info to be passed around
handles.filename = varargin{1};         %filename
handles.event_array = varargin{2};      %event_array
handles.xyzs_id = varargin{3};          %xyzs_id
handles.xyzs_id_columns = varargin{4};  %xyzs_id_columns

image_info = imfinfo(handles.filename);

%Find relevant division information
nume = size(handles.event_array,1); % Number of events
min_frame = 5; % Remove this once we finalize order for determining division correction order

% Initialize vector containing division event numbers
div_event = zeros(nume,1);
for r=1:nume
    if handles.event_array{r}(1,6) == 1 && handles.event_array{r}(1,4) >= min_frame; % If the event is a division, store that event number
        div_event(r) = r;
    end
end

div_event(div_event==0) = [];

handles.div_event = div_event;    

%Check if this is the first time the manual GUI has been called
if size(varargin,2) <=4    
    %Pull out first division event's information
    handles.currentEventIterator = 1;
    handles.currentEventNum = div_event(handles.currentEventIterator);
    handles.currentEventAdvance = 1;
else
   %Set up info from previous run
   handles.currentEventIterator = varargin{5}; 
   handles.currentEventNum = varargin{6};
end

%Set up event counters for user
total_divs = size(handles.div_event,1);
set(handles.event_total_text, 'string', num2str(total_divs))
set(handles.event_text, 'string', num2str(handles.currentEventIterator))

handles.currentEvent.currentEventArray = handles.event_array{handles.currentEventNum,1};
handles.currentEvent.eventSize = size(handles.currentEvent.currentEventArray);
handles.currentEvent.startFrame = handles.currentEvent.currentEventArray(1,4);
handles.currentEvent.endFrame = handles.currentEvent.currentEventArray(handles.currentEvent.eventSize(1),4);

frames_before_div = 5; % # frames before division to plot
frames_after_div = 15; % # frames after division to plot

handles.currentEvent.startFrame = handles.currentEvent.startFrame-frames_before_div;
handles.currentEvent.endFrame = handles.currentEvent.endFrame+frames_after_div;

%Make sure that the frame number is within the video frame number bounds:
if handles.currentEvent.startFrame < 1
    handles.currentEvent.startFrame = 1;
end

if handles.currentEvent.endFrame > size(image_info,1)
    handles.currentEvent.endFrame = size(image_info,1);
end

%Set slider bar info based on first event
sliderMin = handles.currentEvent.startFrame;
sliderMax = handles.currentEvent.endFrame; 

set(handles.video_slider, 'Min', sliderMin);
set(handles.video_slider, 'Max', sliderMax);
set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
set(handles.video_slider, 'Value', sliderMin);

%Set up radio button controls:
% First set default radio button to inivisible button in panel - makes it 
% so that the user input is always recorded
set(handles.division_panel,'SelectedObject', handles.radio3);  
set(handles.division_panel, 'SelectionChangeFcn',  @record_division);

%Set up arrows to advance/reverse events
set(handles.button1, 'FontName', 'Lucida Sans Unicode', 'String', char(60))
set(handles.button2, 'FontName', 'Lucida Sans Unicode', 'String', char(62))

%Set-up video image information:
handles_axis = handles.video_axes;
set(handles.figure1, 'CurrentAxes', handles_axis);
division_plot(handles.currentEvent.startFrame, handles)

% Choose default command line output for division_gui
handles.output = hObject;
if size(varargin,2) <=4
    handles.output = ones(total_divs,1);
else
    handles.output = varargin{7};
end

%Set-up pause button output default
handles.pause_output = 0;

% Update handles structure
guidata(hObject, handles);

%UIWAIT makes division_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = division_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.pause_output;
varargout{3} = handles.div_event;
delete(handles.figure1)


% --- Executes on slider movement.
function video_slider_Callback(hObject, eventdata, handles)
% hObject    handle to video_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get current slider frame number
currentSliderStep = round(get(handles.video_slider, 'Value'));

%Edit plot according to slider position
division_plot(currentSliderStep, handles)

% --- Executes during object creation, after setting all properties.
function video_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to video_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function division_plot(frame_num, handles)

window_buffer = 50; % Pixels buffer for plotting window

image_info = imfinfo(handles.filename);
color_mat = ['b', 'g', 'r', 'y', 'm', 'c']; % Track trajectory colors for video plotting

event_info = handles.event_array{handles.currentEventNum,1};

% Identify dividing cells
unique_cells = nonzeros(unique(event_info(:,2:3)));
num_cells = length(unique_cells);
cell_array = cell(num_cells,1);

% For each cell, find its division event xyzs_id information and identify
% a video window for viewing
for k = 1:num_cells
    cell_pos = handles.xyzs_id(:,handles.xyzs_id_columns) == unique_cells(k) & handles.xyzs_id(:,handles.xyzs_id_columns-1) >= handles.currentEvent.startFrame & handles.xyzs_id(:,handles.xyzs_id_columns-1) <= handles.currentEvent.endFrame;
    cell_mat = handles.xyzs_id(cell_pos,:);
    cell_sort = sortrows(cell_mat,handles.xyzs_id_columns-1);
    cell_array{k} = cell_sort;

    if isempty(cell_sort)
        continue
    end

    % Find minimum and maximum x/y positions for cells involved in the
    % division
    min_xpos(k) = floor(min(cell_sort(:,1)));
    max_xpos(k) = ceil(max(cell_sort(:,1)));
    min_ypos(k) = floor(min(cell_sort(:,2)));
    max_ypos(k) = ceil(max(cell_sort(:,2)));

end

% Identify overall video window for both dividing cells
xmin = max(min(min_xpos)-window_buffer,1);
ymin = max(min(min_ypos)-window_buffer,1);
xmax = min(max(max_xpos)+window_buffer,image_info(1).Width);
ymax = min(max(max_ypos)+window_buffer,image_info(1).Width);

% Plot the full image and incorporate all relevant cell ID values
A = imread(handles.filename, frame_num, 'Info', image_info);
A=imadjust(A);
imagesc(A);
colormap('gray');
title(['Frame ' num2str(frame_num)]);   
hold on

% Change to division event window
axis([xmin xmax ymin ymax])

% Create array for holding cell index information (from xyzs_id)
bool = cell(1,num_cells);

% Loop to find and plot information for all cells involved in
% division event
for j = 1:num_cells
    bool{j} = (cell_array{j}(:,handles.xyzs_id_columns-1)<=frame_num) & (cell_array{j}(:,handles.xyzs_id_columns-1)>=handles.currentEvent.startFrame);
    plot(cell_array{j}(bool{j},1), cell_array{j}(bool{j},2),'marker', 'o', 'color', color_mat(j));
    hold on
end

hold off


function record_division(hObject, eventdata)
handles = guidata(hObject);
switch get(eventdata.NewValue, 'Tag')
    case 'radio1'
        handles.output(handles.currentEventIterator,1) = 1;
        handles.currentEventAdvance = 1;
        
        %Update all of the guidata and run process again
        guidata(hObject,handles)
        
        %Proceed to next event
        change_event(hObject, eventdata, handles)        
    case 'radio2'
        handles.output(handles.currentEventIterator,1) = 0;
        handles.currentEventAdvance = 1;
        
        %Update all of the guidata and run process again
        guidata(hObject,handles)
                
        %Proceed to next event
        change_event(hObject, eventdata, handles)
end
 
function change_event(hObject, eventdata, handles)

%Check which division event the user is currently working on
if handles.currentEventIterator >= size(handles.div_event,1)
    %output structure and close GUI
    uiresume(handles.figure1)
    
%Update event number and set up GUI information    
else
    %Check whether the user has hit a button to move between events
    if handles.currentEventAdvance == 1
        handles.currentEventIterator = handles.currentEventIterator + 1;
        %Update event number that user sees
        set(handles.event_text, 'string', num2str(handles.currentEventIterator))
    elseif handles.currentEventAdvance == 0 && handles.currentEventIterator > 1
        handles.currentEventIterator = handles.currentEventIterator - 1;
        %Update event number that user sees
        set(handles.event_text, 'string', num2str(handles.currentEventIterator))
    end
    
    handles.currentEventNum = handles.div_event(handles.currentEventIterator);
    handles.currentEvent.currentEventArray = handles.event_array{handles.currentEventNum,1};
    handles.currentEvent.eventSize = size(handles.currentEvent.currentEventArray);
    handles.currentEvent.startFrame = handles.currentEvent.currentEventArray(1,4);
    handles.currentEvent.endFrame = handles.currentEvent.currentEventArray(handles.currentEvent.eventSize(1),4);
    
    %Get image information for new event
    image_info = imfinfo(handles.filename);
    frames_before_div = 5; % # frames before division to plot
    frames_after_div = 15; % # frames after division to plot

    handles.currentEvent.startFrame = handles.currentEvent.startFrame-frames_before_div;
    handles.currentEvent.endFrame = handles.currentEvent.endFrame+frames_after_div;

    %Make sure that the frame number is within the video frame number bounds:
    if handles.currentEvent.startFrame < 1
        handles.currentEvent.startFrame = 1;
    end

    if handles.currentEvent.endFrame > size(image_info,1)
        handles.currentEvent.endFrame = size(image_info,1);
    end

    %Set slider bar info based on first event
    sliderMin = handles.currentEvent.startFrame;
    sliderMax = handles.currentEvent.endFrame; 

    set(handles.video_slider, 'Min', sliderMin);
    set(handles.video_slider, 'Max', sliderMax);
    set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
    set(handles.video_slider, 'Value', sliderMin);

    %Set-up video image information:
    division_plot(handles.currentEvent.startFrame, handles)
    
    %Set up radio button controls:
    % First set default radio button to inivisible button in panel - makes it 
    % so that the user input is always recorded
    set(handles.division_panel,'SelectedObject', handles.radio3);  
    set(handles.division_panel, 'SelectionChangeFcn',  @record_division);

    guidata(hObject,handles)

end

% --- Executes on button press in button1.
function button2_Callback(hObject, eventdata, handles)
% hObject    handle to button1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in button2.
handles.currentEventAdvance = 1;

%Note: By default, if the user hits the next arrow without selecting yes or
%no, the event will automatically be declared a division.

%Proceed to next event
change_event(hObject, eventdata, handles)

function button1_Callback(hObject, eventdata, handles)
% hObject    handle to button2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.currentEventAdvance = 0;

%Go to previous event
change_event(hObject, eventdata, handles)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%If the user closes the GUI window without finishing all events, it will 
%still be waiting for user input. Resume Matlab and output results. 
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes during object creation, after setting all properties.
function video_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to video_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate video_axes


% --- Executes on button press in add_frames.
function add_frames_Callback(hObject, eventdata, handles)
% hObject    handle to add_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get slider information from before button push
sliderMin = handles.currentEvent.startFrame;
sliderMax = handles.currentEvent.endFrame; 
currentSliderStep = round(get(handles.video_slider, 'Value'));

%Add an additional 5 frames to the beginning of the video
handles.currentEvent.startFrame = handles.currentEvent.startFrame-5;

%Make sure that the frame number is within the video frame number bounds:
if handles.currentEvent.startFrame < 1
    handles.currentEvent.startFrame = 1;
end

if sliderMin ~= handles.currentEvent.startFrame
    %Update slider bar info
    sliderMin = handles.currentEvent.startFrame;
    set(handles.video_slider, 'Min', sliderMin);
    set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
    set(handles.video_slider, 'Value', currentSliderStep);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in pause_button.
function pause_button_Callback(hObject, eventdata, handles)
% hObject    handle to pause_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Pause button is designed so that the user can continue manual assessment 
%at a later timepoint. Note: if the user decides to pause the manual
%analysis, any user input changes will be recorded and saved but not
%incorporated into the current run. 

%Update the pause status so that the program knows to continue without user
%changes.
handles.pause_output = 1;

%Have the user specify a name for the file and a directory to save it in
[new_filename, pathname] = uiputfile();

%Save the file and close the GUI window
div_output = handles.output;
pause_status = handles.pause_output;
filename = handles.filename;
event_array = handles.event_array;
xyzs_id = handles.xyzs_id;
xyzs_id_columns = handles.xyzs_id_columns; 
currentEventIterator = handles.currentEventIterator;
currentEventNum = handles.currentEventNum;
save([pathname,new_filename], 'div_output', 'pause_status', 'filename', 'event_array', 'xyzs_id', 'xyzs_id_columns', 'currentEventIterator', 'currentEventNum')
uiresume(handles.figure1)
