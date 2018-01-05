function varargout = merging_gui(varargin)
% MERGING_GUI MATLAB code for merging_gui.fig
%
% To run: 
% [merging_output, pause_output] = merging_gui(filename, multi_event_array, xyzs_id, xyzs_id_columns);
%
% To load previous file:
% [merging_output, pause_output] = merging_gui(filename, multi_event_array, xyzs_id, xyzs_id_columns, currentEventNum, merging_output, output_idx);
%
% Megan Brasch, Richard Baker  10/30/13
%
%      MERGING_GUI, by itself, creates a new MERGING_GUI or raises the existing
%      singleton*.
%
%      H = MERGING_GUI returns the handle to a new MERGING_GUI or the handle to
%      the existing singleton*.
%
%      MERGING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MERGING_GUI.M with the given input arguments.
%
%      MERGING_GUI('Property','Value',...) creates a new MERGING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before merging_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to merging_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help merging_gui

% Last Modified by GUIDE v2.5 09-Jan-2014 15:52:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @merging_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @merging_gui_OutputFcn, ...
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


% --- Executes just before merging_gui is made visible.
function merging_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to merging_gui (see VARARGIN)

% Allow relevant info to be passed around
handles.filename = varargin{1};         %filename
handles.eventinfo = varargin{2};        %multi_event_array
handles.xyzs_id = varargin{3};          %xyzs_id
handles.xyzs_id_columns = varargin{4};  %xyzs_id_columns

%Check if currentEventNum exists - if not, set it equal to the first event
if size(varargin,2) <= 4
    handles.currentEventNum = 1;
else
    handles.currentEventNum = varargin{5};
end

%Set up event counter for user
total_events = size(handles.eventinfo,1);
set(handles.event_total_text, 'string', num2str(total_events))
set(handles.event_num_text, 'string', num2str(handles.currentEventNum))

handles.currentEvent.eventnum = handles.eventinfo{handles.currentEventNum,1};
handles.currentEvent.cellvec = handles.eventinfo{handles.currentEventNum,2};
handles.currentEvent.startframe = handles.eventinfo{handles.currentEventNum,3};
handles.currentEvent.endframe = handles.eventinfo{handles.currentEventNum,4}+10; %Add an extra 10 frames to the video

%Set slider bar info based on first event
sliderMin = handles.currentEvent.startframe;
sliderMax = handles.currentEvent.endframe; 

%Check that max does not exceed size of video
image_info = imfinfo(handles.filename);

if sliderMax>size(image_info,1)
    sliderMax = size(image_info,1);
    handles.currentEvent.endframe = size(image_info,1);
end

set(handles.video_slider, 'Min', sliderMin);
set(handles.video_slider, 'Max', sliderMax);
set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
set(handles.video_slider, 'Value', sliderMin);

%Find all cells of interest in first event's window:
xyzs_id = handles.xyzs_id;
ID_vec = handles.currentEvent.cellvec';
xyzs_id_columns = handles.xyzs_id_columns;
bool = zeros(length(xyzs_id),1);
for s=1:length(handles.currentEvent.cellvec)
    bool = bool == 1 | xyzs_id(:,xyzs_id_columns) == ID_vec(s);
end

%Area to look for IDs
xmin2 = min(xyzs_id(bool,1))-5;
ymin2 = min(xyzs_id(bool,2))-5;
xmax2 = max(xyzs_id(bool,1))+5;
ymax2 = max(xyzs_id(bool,2))+5;

%Find all cells present in first and last frame; ignore cells that show up
%but do not appear in either the selected first or last frame. 
bool_vid_cells = (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.endframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2)...
    | (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.startframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2);

full_cell_vec = unique(xyzs_id(bool_vid_cells, xyzs_id_columns));
handles.currentEvent.fullCellVec = full_cell_vec;

%Set-up check for skipping events with more than 24 cells present in a
%complex merging event
handles.no_event = 0;

%Check if extra cell uicontrols need to be visible:
if size(handles.currentEvent.fullCellVec,1)>3
    vis_output = visibility_func(hObject, eventdata, handles);
end

%Set up video plot
ID_vec = handles.currentEvent.fullCellVec';
handles_axis = handles.video_axes;
axes(handles_axis)
[~, ~, ~, ~,  vid_color_mat] = event_video_plot(xyzs_id,xyzs_id_columns, ID_vec, handles.filename, sliderMin, handles_axis, handles);
handles.currentEvent.vid_color_mat = vid_color_mat;

%Set up first frame plot
handles_axis = handles.start_frame_axes;
axes(handles_axis)
[startFrameVec, startFrameNumbers] = event_video_plot(handles.xyzs_id,handles.xyzs_id_columns, ID_vec, handles.filename, sliderMin, handles_axis, handles);
handles.startFrameVec = startFrameVec;
handles.startFrameNumbers = startFrameNumbers;

%Set up end frame plot
handles_axis = handles.end_frame_axes;
axes(handles_axis)
[~, ~, endFrameVec, endFrameLetters] = event_video_plot(handles.xyzs_id,handles.xyzs_id_columns, ID_vec, handles.filename, sliderMax, handles_axis, handles);
handles.endFrameVec = endFrameVec;
handles.endFrameLetters = endFrameLetters;

%Set text under video to display initial correct statement
set(handles.multiBodyText, 'String', 'Multi-Body Interaction', 'ForegroundColor', 'Red')

%Set initial values for edit boxes based on track information
[idMap] = set_values(handles);
handles.currentEvent.idMap = idMap;

%Set-up number values so that they are visible in first frame:
%Figure out which IDs are present in the first frame
bool_cells = handles.currentEvent.idMap(:,3) ~= 0;
cell_ids = handles.currentEvent.idMap(bool_cells,3);

%Pull out color information for plotting
colors = handles.currentEvent.vid_color_mat(bool_cells);

%Set video axes
handles_axis = handles.video_axes;
axes(handles_axis)

for q = 1:size(cell_ids,1)
    %Identify each cell position
    bool_id = (handles.xyzs_id(:,handles.xyzs_id_columns-1) == handles.currentEvent.startframe & handles.xyzs_id(:, handles.xyzs_id_columns) == cell_ids(q));

    %Plot each track with its corresponding number value from first
    %frame still image
    text(handles.xyzs_id(bool_id,1),handles.xyzs_id(bool_id,2),num2str(handles.startFrameNumbers(q)),'FontSize', 14, 'Color', colors(q));  
end

% Choose output for merging_gui
handles.output = hObject;

if size(varargin,2) <= 4
    handles.output = zeros(size(handles.currentEvent.fullCellVec,1),4); %Columns: first frame, last frame, original ID number, updated ID number
    handles.output_idx = 1;     %initialize index used to organize output
else
    handles.output = varargin{6};
    handles.output_idx = varargin{7};
end

%Set-up pause button
handles.pause_output = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes merging_gui wait for user response (see UIRESUME)
uiwait(handles.merging_gui);


% --- Outputs from this function are returned to the command line.
function varargout = merging_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.pause_output;
delete(hObject)


% --- Executes on slider movement.
function video_slider_Callback(hObject, eventdata, handles)
% hObject    handle to video_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get current slider frame number
currentSliderStep = round(get(handles.video_slider, 'Value'));

%Edit plot according to slider position
handles_axis = handles.video_axes;
axes(handles_axis)
ID_vec = handles.currentEvent.fullCellVec';
event_video_plot(handles.xyzs_id,handles.xyzs_id_columns, ID_vec, handles.filename, currentSliderStep, handles_axis, handles);

%Set new_cell_text to be invisible by default
set(handles.new_cell_text,'Visible','Off')

%Check if the user is viewing the first or last frame - show
%numbers/letters accordingly. 

if currentSliderStep == handles.currentEvent.startframe
    
    %Figure out which IDs are present in the first frame
    bool_cells = handles.currentEvent.idMap(:,3) ~= 0;
    cell_ids = handles.currentEvent.idMap(bool_cells,3);
    
    %Pull out color information for plotting
    colors = handles.currentEvent.vid_color_mat(bool_cells);
    
    for q = 1:size(cell_ids,1)
        %Identify each cell position
        bool_id = (handles.xyzs_id(:,handles.xyzs_id_columns-1) == currentSliderStep & handles.xyzs_id(:, handles.xyzs_id_columns) == cell_ids(q));
        
        %Plot each track with its corresponding number value from first
        %frame still image
        text(handles.xyzs_id(bool_id,1),handles.xyzs_id(bool_id,2),num2str(handles.startFrameNumbers(q)),'FontSize', 14, 'Color', colors(q));  
    end
    
elseif currentSliderStep == handles.currentEvent.endframe
    
    %Figure out which IDs are present in the last frame
    bool_cells = handles.currentEvent.idMap(:,4) ~= 0;
    cell_ids = handles.currentEvent.idMap(bool_cells,4);
    
    %Pull out color information for plotting
    colors = handles.currentEvent.vid_color_mat(bool_cells);
    
    for q = 1:size(cell_ids,1)
        %Identify each cell position
        bool_id = (handles.xyzs_id(:,handles.xyzs_id_columns-1) == currentSliderStep & handles.xyzs_id(:, handles.xyzs_id_columns) == cell_ids(q));
        
        %Plot each track with its corresponding letter value from final
        %frame still image
        text(handles.xyzs_id(bool_id,1),handles.xyzs_id(bool_id,2),char(handles.endFrameLetters(q)),'FontSize', 14, 'Color', colors(q));  
    end
end

%Check if a new cell ID has appeared 
if sum(handles.currentEvent.idMap(:,3) == 0) >= 1
    cell_ids = zeros(sum(handles.currentEvent.idMap(:,3) == 0),3);   % cell id, first frame of appearance, xyzs_id reference row
    
    %Figure out which IDs are missing in the first frame
    bool_cells = handles.currentEvent.idMap(:,3) == 0;
    cell_ids(:,1) = handles.currentEvent.idMap(bool_cells,4);
    
    %Find the first frame that each cell is present within the event frame
    %bounds
    for s = 1:sum(handles.currentEvent.idMap(:,3) == 0)
       current_id = cell_ids(s,1); 
       bool_frame = handles.xyzs_id(:,handles.xyzs_id_columns) == current_id & handles.xyzs_id(:,handles.xyzs_id_columns-1) >= handles.currentEvent.startframe & handles.xyzs_id(:,handles.xyzs_id_columns-1) <= handles.currentEvent.endframe;
       xyzs_row_idx = find(bool_frame ~= 0 ,1,'first');
       first_frame = handles.xyzs_id(xyzs_row_idx, handles.xyzs_id_columns-1);
       cell_ids(s,2) = first_frame;
       cell_ids(s,3) = xyzs_row_idx;
    end
    
    %Check if the current slider step = any of the first frame values.
    %Label plot based on assigned number value (col 5 of match_up_mat)
    bool_new_cells = currentSliderStep == cell_ids(:,2);
    if sum(bool_new_cells) > 0 
        frame_ids = cell_ids(bool_new_cells,1);
        for t = 1:size(frame_ids,1)
            cell_id_idx = cell_ids(:,1) == frame_ids(t);
            [row, ~] = find(handles.currentEvent.idMap(:,4) == frame_ids(t));
            hold on
            text(handles.xyzs_id(cell_ids(cell_id_idx,3),1),handles.xyzs_id(cell_ids(cell_id_idx,3),2),num2str(handles.currentEvent.idMap(row,5)),'FontSize', 14, 'Color', handles.currentEvent.vid_color_mat(row));  
        end
        hold off
        set(handles.new_cell_text,'Visible','On')
    end
    
end

%Edit text underneath video window to let user know about state of event
if currentSliderStep <= (handles.currentEvent.endframe - 10)
    set(handles.multiBodyText, 'String', 'Multi-Body Interaction', 'ForegroundColor', 'Red')
else
    set(handles.multiBodyText, 'String', 'Post Multi-Body Interaction', 'ForegroundColor', 'Green')
end

% --- Executes during object creation, after setting all properties.
function video_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to video_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Terminate the manual analysis without evaluating every event - no more
%user input required.
uiresume(handles.merging_gui)

% --- Executes on button press in accept_button.
function accept_button_Callback(hObject, eventdata, handles)
% hObject    handle to accept_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check that the event can be processed (<25 cell present in interaction
%space)
if handles.no_event == 0 
    %Save user input to output structure
    current_event = handles.currentEventNum;
    [new_ids] = create_new_cell_vec(handles);

    %Update output structure to reflect user input
    handles.currentEvent.new_ids = new_ids;
    output_idx = handles.output_idx;
    output_idx = output_idx+size(new_ids,1);

    handles.output(handles.output_idx:output_idx-1, 1) = handles.currentEvent.startframe;
    handles.output(handles.output_idx:output_idx-1, 2) = handles.currentEvent.endframe;
    handles.output(handles.output_idx:output_idx-1, 3) = handles.currentEvent.fullCellVec;
    handles.output(handles.output_idx:output_idx-1, 4) = new_ids;    %Cell IDs based on user input
    handles.output_idx = output_idx;
    
else 
    current_event = handles.currentEventNum;
    handles.no_event = 0;
end 

%Once saved, update the event information
current_event = current_event+1;

if current_event > size(handles.eventinfo,1)
    %resume Matlab - no longer take user input
    uiresume(handles.merging_gui)

%Update event number and set up GUI information    
else
    %Update event counter for user
    set(handles.event_num_text, 'string', num2str(current_event))
    
    %Obtain new event information
    handles.currentEventNum = current_event;
    handles.currentEvent.eventnum = handles.eventinfo{handles.currentEventNum,1};
    handles.currentEvent.cellvec = handles.eventinfo{handles.currentEventNum,2};
    handles.currentEvent.startframe = handles.eventinfo{handles.currentEventNum,3};
    handles.currentEvent.endframe = handles.eventinfo{handles.currentEventNum,4}+10; %Add an extra ten frames to the video

    %Set slider bar info based on first event
    sliderMin = handles.currentEvent.startframe;
    sliderMax = handles.currentEvent.endframe; 
    
    %Check that max does not exceed size of video
    image_info = imfinfo(handles.filename);

    if sliderMax>size(image_info,1)
        sliderMax = size(image_info,1);
        handles.currentEvent.endframe = size(image_info,1);
    end
    
    %Set slider values
    set(handles.video_slider, 'Min', sliderMin);
    set(handles.video_slider, 'Max', sliderMax);
    set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
    set(handles.video_slider, 'Value', sliderMin);
    
    %Find all cells of interest in new event's window:
    xyzs_id = handles.xyzs_id;
    ID_vec = handles.currentEvent.cellvec';
    xyzs_id_columns = handles.xyzs_id_columns;
    bool = zeros(length(xyzs_id),1);
    for s=1:length(handles.currentEvent.cellvec)
        bool = bool == 1 | xyzs_id(:,xyzs_id_columns) == ID_vec(s);
    end

    %Define area to look for IDs
    xmin2 = min(xyzs_id(bool,1))-5;
    ymin2 = min(xyzs_id(bool,2))-5;
    xmax2 = max(xyzs_id(bool,1))+5;
    ymax2 = max(xyzs_id(bool,2))+5;

    %Find all cells present in first and last frame; ignore cells that show
    %up but do not appear in either the selected first or last frame. 
    bool_vid_cells = (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.endframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2)...
    | (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.startframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2);

    full_cell_vec = unique(xyzs_id(bool_vid_cells, xyzs_id_columns));
    handles.currentEvent.fullCellVec = full_cell_vec;

    %Update visibility for text/edit boxes for new event
    vis_output = visibility_func(hObject, eventdata, handles);
    
    %Update handles structure if events have been skipped
    if vis_output ~= 0
        
        current_event = vis_output;
    
        %Update event counter for user
        set(handles.event_num_text, 'string', num2str(current_event))

        %Obtain new event information
        handles.currentEventNum = current_event;
        handles.currentEvent.eventnum = handles.eventinfo{handles.currentEventNum,1};
        handles.currentEvent.cellvec = handles.eventinfo{handles.currentEventNum,2};
        handles.currentEvent.startframe = handles.eventinfo{handles.currentEventNum,3};
        handles.currentEvent.endframe = handles.eventinfo{handles.currentEventNum,4}+10; %Add an extra ten frames to the video

        %Set slider bar info based on first event
        sliderMin = handles.currentEvent.startframe;
        sliderMax = handles.currentEvent.endframe; 

        %Check that max does not exceed size of video
        image_info = imfinfo(handles.filename);

        if sliderMax>size(image_info,1)
            sliderMax = size(image_info,1);
            handles.currentEvent.endframe = size(image_info,1);
        end

        %Set slider values
        set(handles.video_slider, 'Min', sliderMin);
        set(handles.video_slider, 'Max', sliderMax);
        set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
        set(handles.video_slider, 'Value', sliderMin);

        %Find all cells of interest in new event's window:
        xyzs_id = handles.xyzs_id;
        ID_vec = handles.currentEvent.cellvec';
        xyzs_id_columns = handles.xyzs_id_columns;
        bool = zeros(length(xyzs_id),1);
        for s=1:length(handles.currentEvent.cellvec)
            bool = bool == 1 | xyzs_id(:,xyzs_id_columns) == ID_vec(s);
        end

        %Define area to look for IDs
        xmin2 = min(xyzs_id(bool,1))-5;
        ymin2 = min(xyzs_id(bool,2))-5;
        xmax2 = max(xyzs_id(bool,1))+5;
        ymax2 = max(xyzs_id(bool,2))+5;

        %Find all cells present in first and last frame; ignore cells that show
        %up but do not appear in either the selected first or last frame. 
        bool_vid_cells = (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.endframe ...
        & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2)...
        | (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.startframe ...
        & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2);

        full_cell_vec = unique(xyzs_id(bool_vid_cells, xyzs_id_columns));
        handles.currentEvent.fullCellVec = full_cell_vec;

        %Update visibility for text/edit boxes for new event
        vis_output = visibility_func(hObject, eventdata, handles);
    
    end
        
    %Set up video plot
    ID_vec = handles.currentEvent.fullCellVec';
    handles_axis = handles.video_axes;
    axes(handles_axis)
    [~, ~, ~, ~,  vid_color_mat] = event_video_plot(xyzs_id,xyzs_id_columns, ID_vec, handles.filename, sliderMin, handles_axis, handles);
    handles.currentEvent.vid_color_mat = vid_color_mat;
    
    %Set up first frame plot
    handles_axis = handles.start_frame_axes;
    axes(handles_axis)
    [startFrameVec, startFrameNumbers] = event_video_plot(handles.xyzs_id,handles.xyzs_id_columns, ID_vec, handles.filename, sliderMin, handles_axis, handles);
    handles.startFrameVec = startFrameVec;
    handles.startFrameNumbers = startFrameNumbers;

    %Set up end frame plot
    handles_axis = handles.end_frame_axes;
    axes(handles_axis)
    [~, ~, endFrameVec, endFrameLetters] = event_video_plot(handles.xyzs_id,handles.xyzs_id_columns, ID_vec, handles.filename, sliderMax, handles_axis, handles);
    handles.endFrameVec = endFrameVec;
    handles.endFrameLetters = endFrameLetters;
 
    %Set text under video to display initial correct statement
    set(handles.multiBodyText, 'String', 'Multi-Body Interaction', 'ForegroundColor', 'Red')
    
    %Set initial values for edit boxes based on track information
    [idMap] = set_values(handles);
    handles.currentEvent.idMap = idMap;
    
    %Set-up number values so that they are visible in first frame:
    %Figure out which IDs are present in the first frame
    bool_cells = handles.currentEvent.idMap(:,3) ~= 0;
    cell_ids = handles.currentEvent.idMap(bool_cells,3);

    %Pull out color information for plotting
    colors = handles.currentEvent.vid_color_mat(bool_cells);
    
    %Set video axes
    handles_axis = handles.video_axes;
    axes(handles_axis)
    
    for q = 1:size(cell_ids,1)
        %Identify each cell position
        bool_id = (handles.xyzs_id(:,handles.xyzs_id_columns-1) == handles.currentEvent.startframe & handles.xyzs_id(:, handles.xyzs_id_columns) == cell_ids(q));

        %Plot each track with its corresponding number value from first
        %frame still image
        text(handles.xyzs_id(bool_id,1),handles.xyzs_id(bool_id,2),num2str(handles.startFrameNumbers(q)),'FontSize', 14, 'Color', colors(q));  
    end
    
    %Update all of the guidata and run process again
    guidata(hObject,handles)
end

%Organize new cell vector IDs based on user inputs
function [new_ids] = create_new_cell_vec(handles)
    %Pull out ID map information from set_values function
    idMap = zeros(size(handles.currentEvent.idMap,1),6);
    idMap(:,1:5) = handles.currentEvent.idMap;
    
    %Obtain the current user inputs and associate them with original
    %matches
    for w = 1:size(idMap,1)
        editbox = ['edit', num2str(w)];
        current_value = get(handles.(editbox), 'String');
        bool_num = idMap(:,1) == w;
        if sum(bool_num) > 0 
            idMap(bool_num,6) = current_value;
        else
            bool_num2 = idMap(:,5) == w;
            idMap(bool_num2,6) = current_value;
        end
    end
    
    new_vec = idMap(:,6); 
    original_vec = idMap(:,2); 
    [bool_id_match,loc] = ismember(new_vec,original_vec); %Note: 48 in new_vec denotes char code for value 0 - when using ismember, this will appear as a zero in the bool_match_up variable
    
    new_ids = zeros(size(new_vec,1),1);
    
    %This can probably be done without a for loop, but easiest way - small
    %matrices (no larger than 24), so it should not effect speed too much
for b = 1:size(new_vec,1)
    %Find IDs that existed in start frame
    if bool_id_match(b) ~= 0 && idMap(loc(b),3) ~= 0
        new_ids(b) = idMap(loc(b),3);
    %Find IDs that did not exist in start frame
    elseif bool_id_match(b) ~= 0 && idMap(loc(b),3) == 0
        new_ids(b) = idMap(loc(b),4);
    %Zero denotes cells that have disappeared during video time frame
    else
        new_ids(b) = 0;
    end
end

%Custom function to plot videos with specific window on different axes
function [startFrameVec, startFrameNumbers, endFrameVec, endFrameLetters,  vid_color_mat] = event_video_plot(xyzs_id,xyzs_id_columns, ID_vec, stackname, frame_num, handles_axis, handles)

frame_padding = 20;

%Set up outputs for video axis calls
vid_color_mat = [];
if handles_axis == handles.video_axes
    startFrameVec = [];
    startFrameNumbers = [];
    endFrameVec = [];
    endFrameLetters =  [];
end

%Identify cells present in ID_vec
bool = zeros(length(xyzs_id),1);
for s=1:length(ID_vec)
    bool = bool == 1 | xyzs_id(:,xyzs_id_columns) == ID_vec(s);
end

%Plot window bounds
xmin = min(xyzs_id(bool,1))-frame_padding;
ymin = min(xyzs_id(bool,2))-frame_padding;
xmax = max(xyzs_id(bool,1))+frame_padding;
ymax = max(xyzs_id(bool,2))+frame_padding;

%Area to look for IDs
xmin2 = min(xyzs_id(bool,1))-5;
ymin2 = min(xyzs_id(bool,2))-5;
xmax2 = max(xyzs_id(bool,1))+5;
ymax2 = max(xyzs_id(bool,2))+5;

%Read in image information and plot on correct axis
image_info = imfinfo(stackname);

A = imread(stackname, frame_num, 'Info', image_info);
A = imadjust(A);
imagesc(A);
colormap gray
hold on     
axis([xmin xmax ymin ymax]);
title(['Frame ' num2str(frame_num)]);   

%Find all cells for a frame
bool_cells = xyzs_id(:,xyzs_id_columns-1) == frame_num & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2;
new_idx = find(bool_cells~=0&bool==1);

%Set up plotting matrices
base_color_mat = ['b', 'g', 'r', 'y', 'm', 'c']; % Track trajectory colors for video plotting
alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];

%All cells may not be present in first/last frame - process each case
%individually. 
if handles_axis == handles.video_axes
    num_colors = size(handles.currentEvent.fullCellVec,1);
    if num_colors <= size(base_color_mat,2)
        vid_color_mat = base_color_mat;
    else
        vid_color_mat = [];
        multiple = floor(num_colors/size(base_color_mat,2));
        remain = rem(num_colors,size(base_color_mat,2));
        for i = 1:multiple
            vid_color_mat = [vid_color_mat,base_color_mat];
        end
        vid_color_mat = [vid_color_mat,base_color_mat(1:remain)];
    end
    handles.currentEvent.vid_color_mat = vid_color_mat;
end

%Find info for start or end frames
if handles_axis == handles.start_frame_axes || handles_axis == handles.end_frame_axes
    full_id_list = zeros(size(new_idx,1), 4);
    full_id_list(:,1) = xyzs_id(new_idx,xyzs_id_columns);    %list of IDs involved in event
    full_id_list(:,2) = xyzs_id(new_idx,1);  %x position
    full_id_list(:,3) = xyzs_id(new_idx,2);  %y position
    full_id_list(:,4) = (1:size(full_id_list,1))';    % numbers for start frame labeling
    full_id_list(:,5) = alphabet(1:size(full_id_list,1))';  % letters for end frame labeling
end

%Record cells present in the starting and ending frames
if handles_axis == handles.start_frame_axes
    startFrameVec = full_id_list(:,1);
    startFrameNumbers = full_id_list(:,4);
    endFrameVec = [];
    endFrameLetters =  [];
elseif handles_axis == handles.end_frame_axes
    endFrameVec = full_id_list(:,1);
    endFrameLetters =  full_id_list(:,5);
    startFrameVec = [];
    startFrameNumbers = [];
end

if handles_axis == handles.video_axes
    %Plot each cell's trajectory as its own color for every frame it is
    %present during this complex interaction
    num_p_frame = size(handles.currentEvent.fullCellVec,1);
    for k=1:num_p_frame
        bool_single_cell = xyzs_id(:,xyzs_id_columns-1) <= frame_num & xyzs_id(:,xyzs_id_columns-1) >= handles.currentEvent.startframe & xyzs_id(:,xyzs_id_columns) == handles.currentEvent.fullCellVec(k);
        plot(xyzs_id(bool_single_cell, 1), xyzs_id(bool_single_cell,2), ['x-', vid_color_mat(k)]);
    end
    hold off
elseif handles_axis == handles.start_frame_axes
    text(full_id_list(:,2),full_id_list(:,3),num2str(full_id_list(:,4)),'FontSize', 14, 'Color', 'y');  
    hold off
else
    text(full_id_list(:,2),full_id_list(:,3),char(full_id_list(:,5)),'FontSize', 14, 'Color', 'y');  
    hold off
end

%Function to customize text and edit box visibility
function vis_output = visibility_func(hObject, eventdata, handles)
    
%Identify total number of cells and total potential number of edit/text boxes    
num_cells = size(handles.currentEvent.fullCellVec,1);
max_cells = 24; %Pre-determined number of maximum cells per interaction event

%Skip event if there are more than 24 cells present in the window - too
%complex to process manually. 
if num_cells <= max_cells
    
    handles.no_event = 0;
    
    %Loop starts at 4 because cells 1-3 will always be present for a
    %multi-body interaction
    for i = 4:num_cells
        textbox = ['text',num2str(i)];
        set(handles.(textbox), 'Visible', 'On')
        editbox = ['edit', num2str(i)];
        set(handles.(editbox), 'Visible', 'On')
    end

    %Turn off all other box and text visibility:
    for j = num_cells+1:max_cells
        textbox = ['text',num2str(j)];
        set(handles.(textbox), 'Visible', 'Off')
        editbox = ['edit', num2str(j)];
        set(handles.(editbox), 'Visible', 'Off')
    end
    
    vis_output = 0; 
else
    %If more than 24 cells are present, too complex to manually analyze.
    %Skip to next event
    handles.no_event = 1;
    
    [event_num] = skip_event(hObject, eventdata, handles);
    
    vis_output = event_num;
    
end

function [idMap] = set_values(handles)
    %Figure out how many cells are present in the start and end frames
    startFrameVec = handles.startFrameVec;
    endFrameVec = handles.endFrameVec;
    endFrameLetters = handles.endFrameLetters;
    fullCellVec = handles.currentEvent.fullCellVec;
    startFrameNumbers = handles.startFrameNumbers;
    
    %Provide an initial 'guess' for edit boxes based on tracking results
    [bool_id,sIdx] = ismember(fullCellVec,startFrameVec);
    [bool_id2,eIdx] = ismember(fullCellVec,endFrameVec);
    eIdx(bool_id2) = endFrameLetters;
    sIdx(bool_id) = startFrameNumbers;

    %Determine what cells are present in the starting vector and the ending
    %cell vectors
    match_up_mat = [sIdx, eIdx];    %col 1 = cells present in startFrameVec, col2 = cells present in endFrameVec
    zero_loc = find(match_up_mat(:,1) == 0);
    
    %Add ID's into matrix for mapping purposes
    sIdx(bool_id) = startFrameVec;
    eIdx(bool_id2) = endFrameVec;
    match_up_mat = [match_up_mat, sIdx, eIdx]; %col 3 = cell IDs from startFrameVec, col 4 = cell IDs from endFrameVec
    match_up_mat(:,5) = 0;  %Reserve a new column for number values of cells missing in first frame 

    %If cells are missing in the first frame, update column five with new cell number 
    if ~isempty(zero_loc)
        for i = 1:size(zero_loc,1)
            if i == 1
                %Establish the new number for cell missing in first frame
                new_num = max(match_up_mat(:,1))+1;
                match_up_mat(zero_loc(i),5) = new_num;
            else
                %If more than one cell is missing in first frame, increase
                %this value
                new_num = new_num+1;
                match_up_mat(zero_loc(i),5) = new_num;
            end
        end
    end
    
    %Update cell Ids
    num_cells = size(fullCellVec,1);
    for i = 1:num_cells
        %Check if value is missing (cell appears after 1st frame)
        zero_chk = match_up_mat(i,1) == 0;
        if zero_chk > 0 
            %Look at the 5th column to process IDs that appear after
            %the starting frame
            cell_number = match_up_mat(i,5);
            editbox = ['edit', num2str(cell_number)];

            if match_up_mat(i,2) ~= 0
                set(handles.(editbox), 'String', char(match_up_mat(i,2)))
            else 
                set(handles.(editbox), 'String', match_up_mat(i,2))
            end
        else 
            %Otherwise, use the number reference in column 1 of the
            %match_up_mat to determine edit box to set
            cell_number = match_up_mat(i,1);
            editbox = ['edit', num2str(cell_number)];

            if match_up_mat(i,2) ~= 0
                set(handles.(editbox), 'String', char(match_up_mat(i,2)))
            else 
                set(handles.(editbox), 'String', match_up_mat(i,2))
            end
        end
    end

    %Print final matrix to be saved in handles structure for additional
    %manipulation
    idMap = match_up_mat;
    
function [event_num] = skip_event(hObject, eventdata, handles)

event_test = 0;
current_event = handles.currentEventNum;

while event_test == 0 && current_event < size(handles.eventinfo,1)

    %Update event counter for user
    set(handles.event_num_text, 'string', num2str(current_event))

    %Obtain new event information
    handles.currentEventNum = current_event;
    handles.currentEvent.eventnum = handles.eventinfo{handles.currentEventNum,1};
    handles.currentEvent.cellvec = handles.eventinfo{handles.currentEventNum,2};
    handles.currentEvent.startframe = handles.eventinfo{handles.currentEventNum,3};
    handles.currentEvent.endframe = handles.eventinfo{handles.currentEventNum,4}+10; %Add an extra ten frames to the video

    %Set slider bar info b7ased on first event
    sliderMin = handles.currentEvent.startframe;
    sliderMax = handles.currentEvent.endframe; 

    %Check that max does not exceed size of video
    image_info = imfinfo(handles.filename);

    if sliderMax>size(image_info,1)
        sliderMax = size(image_info,1);
        handles.currentEvent.endframe = size(image_info,1);
    end

    %Set slider values
    set(handles.video_slider, 'Min', sliderMin);
    set(handles.video_slider, 'Max', sliderMax);
    set(handles.video_slider, 'SliderStep', [1/(sliderMax-sliderMin) 0.25]);
    set(handles.video_slider, 'Value', sliderMin);

    %Find all cells of interest in new event's window:
    xyzs_id = handles.xyzs_id;
    ID_vec = handles.currentEvent.cellvec';
    xyzs_id_columns = handles.xyzs_id_columns;
    bool = zeros(length(xyzs_id),1);
    for s=1:length(handles.currentEvent.cellvec)
        bool = bool == 1 | xyzs_id(:,xyzs_id_columns) == ID_vec(s);
    end

    %Define area to look for IDs
    xmin2 = min(xyzs_id(bool,1))-5;
    ymin2 = min(xyzs_id(bool,2))-5;
    xmax2 = max(xyzs_id(bool,1))+5;
    ymax2 = max(xyzs_id(bool,2))+5;

    %Find all cells present in first and last frame; ignore cells that show
    %up but do not appear in either the selected first or last frame. 
    bool_vid_cells = (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.endframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2)...
    | (xyzs_id(:,xyzs_id_columns-1) == handles.currentEvent.startframe ...
    & xyzs_id(:,1) > xmin2 & xyzs_id(:,1) < xmax2 & xyzs_id(:,2) > ymin2 & xyzs_id(:,2) < ymax2);

    full_cell_vec = unique(xyzs_id(bool_vid_cells, xyzs_id_columns));
    handles.currentEvent.fullCellVec = full_cell_vec;

    if size(handles.currentEvent.fullCellVec,1) > 24
       current_event = current_event+1;
    else
        event_test = 1;
        
        %Identify new event number 
        event_num = current_event;
    end
    
end

    % --- Executes during object creation, after setting all properties.
function start_frame_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_frame_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function end_frame_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_frame_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close merging_gui.
function merging_gui_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to merging_gui (see GCBO)
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


% --- Executes on button press in pause_button.
function pause_button_Callback(hObject, eventdata, handles)
% hObject    handle to pause_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Update the pause status so that the program knows to continue without user
%changes.
handles.pause_output = 1;

%Have the user specify a name for the file and a directory to save it in
[new_filename, pathname] = uiputfile();

%Save the file and close the GUI window
merging_output = handles.output;
pause_status = handles.pause_output;
filename = handles.filename;
multi_event_array = handles.eventinfo;
xyzs_id = handles.xyzs_id;
xyzs_id_columns = handles.xyzs_id_columns; 
currentEventNum = handles.currentEventNum;
output_idx = handles.output_idx;
save([pathname,new_filename], 'merging_output', 'pause_status', 'filename', 'multi_event_array', 'xyzs_id', 'xyzs_id_columns','currentEventNum', 'output_idx')
uiresume(handles.merging_gui)


%% Code for all of the following edit boxes is identical:
% 24 edit boxes were included as a default. Their visibility is toggled
% on/off depending on the number of cells present in the first/last frame

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

currentVal = get(hObject,'String');
%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit1, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit2, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit3, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit4, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit5, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit6, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit7, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit8, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit9, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit10, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit11, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit12, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit13, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit14, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit15, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit16, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit17, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit18, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit19, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit20, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit21, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit22, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit23, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');


function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVal = get(hObject,'String');

%Check that the user entered a proper value
if (ischar(currentVal) && ismember(currentVal,handles.endFrameLetters)) || str2double(currentVal) == 0
    set(handles.edit24, 'String', currentVal);
    set(handles.warning_text, 'Visible', 'off')
else
    set(handles.warning_text, 'String', 'You must input a capital letter present in the final frame or a value of 0 for a cell that disappears.', 'Visible', 'on')
end

% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Make this edit box invisible unless there are more than 3 cells involved
%in this complex merging event
set(hObject,'Visible', 'off');

