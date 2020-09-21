function varargout = HotSpot(varargin)
% HOTSPOT MATLAB code for HotSpot.fig
%      HOTSPOT, by itself, creates a new HOTSPOT or raises the existing
%      singleton*.
%
%      H = HOTSPOT returns the handle to a new HOTSPOT or the handle to
%      the existing singleton*.
%
%      HOTSPOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HOTSPOT.M with the given input arguments.
%
%      HOTSPOT('Property','Value',...) creates a new HOTSPOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HotSpot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HotSpot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HotSpot

% Last Modified by GUIDE v2.5 21-Sep-2020 15:55:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HotSpot_OpeningFcn, ...
                   'gui_OutputFcn',  @HotSpot_OutputFcn, ...
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


% --- Executes just before HotSpot is made visible.
function HotSpot_OpeningFcn(hObject, eventdata, handles, varargin)

%- Overlay the different analysis tools
pos_hs = get(handles.panel_hotspot,'Position');
set(handles.panel_coloc,'Position',pos_hs);

%- Plot status
handles.status_zoom = 0;
handles.status_pan  = 0;
    
handles.status_first_plot = 1;

%- File-identifier to open rna images
handles.file_ident.status = 0;
handles.file_ident.protein = 'w1GFP';
handles.file_ident.rna     = 'w2Cy3';
    
%- Dummy FQ object
handles.img_rna     = FQ_img;
handles.img_prot    = FQ_img;
handles.spot_counter = 1;
handles = translation_init(hObject, eventdata, handles);

%- Change size of averaging
handles.img_prot.settings.avg_spots.crop.xy = 3;
handles.img_prot.settings.avg_spots.crop.z  = 1;

%- Link-axes
linkaxes([handles.axes_protein,handles.axes_rna],'xy');

%- Default analysis is HotSpots
popup_analysis_Callback(hObject, eventdata, handles);


% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = HotSpot_OutputFcn(hObject, eventdata, handles) 
%varargout{1} = handles.output;


function handles = translation_init(hObject, eventdata, handles);

%- Reinit for coloc
handles.cell_counter = 1;
handles.cell_prop    = struct('name','','x',[],'y',[],'x_prot',[],'y_prot',[],'x_rna',[],'y_rna',[]);

set(handles.listbox_cl_cells,'String','');
set(handles.listbox_cl_cells,'Value',1);

%- Reinit for hotspots
handles.spot_counter = 1;
handles.translation = struct('name','','x',[],'y',[],'z',[],'x_prot',[],'y_prot',[],'z_prot',[]);
set(handles.listbox_translation,'String','');
set(handles.listbox_translation,'Value',1);

handles.x_prot = [];
handles.y_prot = [];
handles.z_prot = [];

handles.x_prot_all = [];
handles.y_prot_all = [];
handles.z_prot_all = [];
handles.ind_all    = [];

guidata(hObject, handles);


%% ==== Define experimental parameters
function menu_settings_Callback(hObject, eventdata, handles)

%- Define input dialog
dlgTitle = 'Analysis parameters';

prompt(1) = {'Crop-size XY [pix]'};
prompt(2) = {'Crop-size Z [slices, 0 for 2D]'};

defaultValue{1} = num2str(handles.img_prot.settings.avg_spots.crop.xy);
defaultValue{2} = num2str(handles.img_prot.settings.avg_spots.crop.z);

userValue = inputdlg(prompt,dlgTitle,1,defaultValue);

if( ~ isempty(userValue))
    handles.img_prot.settings.avg_spots.crop.xy = str2double(userValue{1});
    handles.img_prot.settings.avg_spots.crop.z  = str2double(userValue{2});   
    guidata(hObject, handles);
end


%==========================================================================
% LOAD IMAGE DATA
%==========================================================================

%= Load protein data
function button_load_protein_Callback(hObject, eventdata, handles)

%- First time calling - define identifiers
if  ~handles.file_ident.status
    
    prompt = {'Protein','RNA'};
    dlg_title = 'Unique channel identifiers';
    num_lines = 1;
    def = {handles.file_ident.protein,handles.file_ident.rna};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        handles.file_ident.protein = answer{1};
        handles.file_ident.rna = answer{2};         
    else
        return
    end
    
end

%- Open protein data
handles.img_prot.reinit;
status_file = handles.img_prot.load_img([],'raw');

if ~status_file; return; end

%- Perform z-projection
handles.img_prot.project_Z('raw','max');

%- Controls for z-stacks
set(handles.text_frame,'String','1');
set(handles.slider_slice,'Value',0)

%- Contrast
handles.img_prot_min  =  double(min(handles.img_prot.raw(:))); 
handles.img_prot_max  =  double(max(handles.img_prot.raw(:))); 
handles.img_prot_diff =  double(handles.img_prot_max-handles.img_prot_min);
handles.img_prot_NZ   =  size(handles.img_prot.raw,3); 

%- Contrast
set(handles.text_contr_prot_min,'String',num2str(round(handles.img_prot_min)));
set(handles.text_contr_prot_max,'String',num2str(round(handles.img_prot_max)));

set(handles.slider_contrast_prot_min,'Value',0)
set(handles.slider_contrast_prot_max,'Value',1)


%- Open rna image
if ~isempty(handles.file_ident.rna)
    default_name = strrep(handles.img_prot.file_names.raw, handles.file_ident.protein, handles.file_ident.rna);  
    
    %- Can name of protein image be obtained by string replacement
    %  Can be the same name if no replacement was performed ... 
    if not(strcmp( handles.file_ident.rna,handles.file_ident.protein)) && ...
       isequal(handles.img_prot.file_names.raw,default_name)
        errordlg('Could not find correct name for rna data by string replacement. see command line for details')
        fprintf('File-name    : %s\n',handles.img_prot.file_names.raw)
        fprintf('Ident-protein: %s\n',handles.file_ident.protein)
        fprintf('Ident-protein: %s\n',handles.file_ident.rna)
        
    else
        
        %- Does protein image file exist?
        if ~exist(fullfile(handles.img_prot.path_names.img ,default_name))
            errordlg('rna image does not exist. see command line for details')
            fprintf('File-name    : %s\n',default_name)
            fprintf('Ident-protein: %s\n',handles.img_prot.path_names.img)
        else
            status_file = handles.img_rna.load_img(fullfile(handles.img_prot.path_names.img ,default_name),'raw');
            
            if status_file
                
                %- Perform z-projection
                handles.img_rna.project_Z('raw','max');

                handles.img_rna_min  =  double(min(handles.img_rna.raw(:))); 
                handles.img_rna_max  =  double(max(handles.img_rna.raw(:))); 
                handles.img_rna_diff =  double(handles.img_rna_max-handles.img_rna_min);
                handles.img_rna_NZ   =  size(handles.img_rna.raw,3); 

                %- Contrasts
                set(handles.text_contr_rna_min,'String',num2str(round(handles.img_rna_min)));
                set(handles.text_contr_rna_max,'String',num2str(round(handles.img_rna_max)));
                
                set(handles.slider_contrast_rna_min,'Value',0)
                set(handles.slider_contrast_rna_max,'Value',1)
                
                handles.file_ident.status = 1;
                
                
            else
               	errordlg('rna image could not be loaded. see command line for details')
                fprintf('File-name    : %s\n',default_name)
                fprintf('Ident-protein: %s\n',handles.img_prot.path_names.img)
            end
        end
    end
end

%- Plot image    
handles = translation_init(hObject, eventdata, handles);
handles.status_first_plot = 1;
plot_image(hObject, eventdata, handles);

%- Save data
guidata(hObject, handles);


%==========================================================================
% Choose analysis
%==========================================================================

function popup_analysis_Callback(hObject, eventdata, handles)
if get(handles.popup_analysis,'Value') == 1
    set(handles.panel_hotspot,'Visible','on')
    set(handles.panel_coloc,'Visible','off')
    
    set(handles.menu_hotspot,'enable','on')
    set(handles.menu_coloc,'enable','off')
    
    set(handles.panel_zstack,'Visible','on')
    
elseif get(handles.popup_analysis,'Value') == 2
    set(handles.panel_hotspot,'Visible','off')
    set(handles.panel_coloc,'Visible','on')
    
    set(handles.menu_hotspot,'enable','off')
    set(handles.menu_coloc,'enable','on')
    
    set(handles.panel_zstack,'Visible','off')
end

plot_image(hObject, eventdata, handles);

%==========================================================================
% Analysis - hotspots
%==========================================================================

%= Define experimental parameters
function menu_hotspot_Callback(hObject, eventdata, handles)


%= Define translation-spot
function button_define_translation_Callback(hObject, eventdata, handles)

[x,y] = ginputax(handles.axes_protein,1);
x = round(x); y=round(y); 
z = str2double(get(handles.text_frame,'String'));

%- Correct position
crop = handles.img_prot.settings.avg_spots.crop;

xmin = (x - crop.xy); 
xmax = (x + crop.xy);

ymin = (y - crop.xy);
ymax = (y + crop.xy);

zmin = (z - crop.z);
zmax = (z + crop.z);
    
%- Cropped image
try
    img_crop = handles.img_prot.raw(ymin:ymax,xmin:xmax,zmin:zmax);
catch 
    errordlg('Spot close to edge of image (in X,Y, or Z)','HotSpot')
    return
end
    
    
%- Find maximum position in this image
[dum, Ixyz] = max(img_crop(:));
[Iy,Ix,Iz] = ind2sub(size(img_crop),Ixyz);

%- Correct original estimates by this offset
x_new = x + (Ix-crop.xy-1);
y_new = y + (Iy-crop.xy-1);
z_new = z + (Iz-crop.z-1);

% x=x_new;y=y_new;z=z_new;
%- Define values of listbox
str_list = get(handles.listbox_translation,'String');
N_trans   = numel(str_list);
ind_trans = N_trans+1;

str_cell = ['TranslationSpot_', num2str(handles.spot_counter)];
str_list{ind_trans} = str_cell;
set(handles.listbox_translation,'String',str_list);        
set(handles.listbox_translation,'Value',ind_trans);

handles.translation(ind_trans).name = str_cell;
handles.translation(ind_trans).x    = x_new;
handles.translation(ind_trans).y    = y_new;
handles.translation(ind_trans).z    = z_new;

%- Save data & plot
handles.spot_counter = handles.spot_counter+1;
guidata(hObject, handles);
listbox_translation_Callback(hObject, eventdata, handles)


%= Delete translation spot
function button_delete_translation_Callback(hObject, eventdata, handles)

    
%- Extract index of highlighted cell
str_list = get(handles.listbox_translation,'String');
ind_sel  = get(handles.listbox_translation,'Value');

%- Delete highlighted cell
str_list(ind_sel) = [];

handles.translation(ind_sel) = [];   
set(handles.listbox_translation,'String',str_list)

%- Make sure that pointer is not outside of defined translation spots
N_str = length(str_list);
if ind_sel > N_str     
    set(handles.listbox_translation,'Value',N_str)
end

%- Save data & plot
guidata(hObject, handles);
listbox_translation_Callback(hObject, eventdata, handles)


%= Delete all translation spot
function button_delete_translation_all_Callback(hObject, eventdata, handles)

%- Ask user to confirm choice
choice = questdlg('Do you really want to delete ALL translation spots?', mfilename, 'Yes','No','No');

if strcmp(choice,'Yes')
    
    %- Save data & plot
    handles = translation_init(hObject, eventdata, handles);
    guidata(hObject, handles);
    listbox_translation_Callback(hObject, eventdata, handles)
    
end
    
%==========================================================================
% Analysis - proteins
%==========================================================================

%= Define individual proteins
function button_define_proteins_Callback(hObject, eventdata, handles)

%- Check if translation spot is defined
ind_trans = get(handles.listbox_translation,'Value');
if isempty(handles.translation(ind_trans).x)
    errordlg('No translation spot defined!','HotSpot')
    return
end

%- Define proteins
x = [];
y = [];

while 1
    [x_loop,y_loop,button,axn] = ginputax(handles.axes_protein,1);
   
    %- Empty position - button pressed
    if isempty(x_loop); break; end
    
    %- Check if in actual axes (and not in axes for RNA)
    if isempty(axn); break; end
    
    %- Make sure not outside of the axes but in the GUI
    if  (x_loop < 0.5  ||  x_loop > (handles.img_prot.dim.X+0.5) && ...
       y_loop < 0.5  ||  y_loop > (handles.img_prot.dim.Y+0.5))
        
        break
    end
    
    %- Save positions
    x=[x;x_loop]; y=[y;y_loop];

    %- Show image
    hold on
        plot(handles.axes_protein,x_loop,y_loop,'g+')
    hold off
 
end

%- Return if no position was defined
if isempty(x); return; end

%- Analyze positions
x = round(x);
y = round(y); 
z = repmat(str2double(get(handles.text_frame,'String')),numel(x),1);

%- Go over spots and find maximum
crop = handles.img_prot.settings.avg_spots.crop;

for iS = 1:numel(x)
    
    %- Cropping window
    xmin = (x(iS) - crop.xy); 
    xmax = (x(iS) + crop.xy);
    
    ymin = (y(iS) - crop.xy);
    ymax = (y(iS) + crop.xy);
    
    zmin = (z(iS) - crop.z);
    zmax = (z(iS) + crop.z);
    
    %- Cropped image
    try
        img_crop = handles.img_prot.raw(ymin:ymax,xmin:xmax,zmin:zmax);
    catch
        listbox_translation_Callback(hObject, eventdata, handles)
        errordlg('Spot close to edge of image (in X,Y, or Z)','HotSpot')
        return
    end
    
    %- Find maximum position in this image
    [dum, Ixyz] = max(img_crop(:));
    [Iy,Ix,Iz]  = ind2sub(size(img_crop),Ixyz);
    
    %- Correct original estimates by this offset
    x_new(iS,1) = x(iS) + (Ix-crop.xy-1);
    y_new(iS,1) = y(iS) + (Iy-crop.xy-1);
    z_new(iS,1) = z(iS) - (Iz-crop.z-1);
    
end

if isempty(handles.translation(ind_trans).x_prot)
    handles.translation(ind_trans).x_prot = round(x_new);
    handles.translation(ind_trans).y_prot = round(y_new);
    handles.translation(ind_trans).z_prot = round(z_new);
    
else
    handles.translation(ind_trans).x_prot = [handles.translation(ind_trans).x_prot;round(x_new)];
    handles.translation(ind_trans).y_prot = [handles.translation(ind_trans).y_prot;round(y_new)];
    handles.translation(ind_trans).z_prot = [handles.translation(ind_trans).z_prot;round(z_new)];
end

%- Add to list
handles.x_prot_all = [handles.x_prot_all;x_new];
handles.y_prot_all = [handles.y_prot_all;y_new];
handles.z_prot_all = [handles.z_prot_all;z_new];
handles.ind_all    = [handles.ind_all;repmat(ind_trans,size(x_new))];

%- Save positions for later
handles.x_prot = x_new;
handles.y_prot = y_new;
handles.z_prot = z_new;

%- Save data & plot
guidata(hObject, handles);
listbox_translation_Callback(hObject, eventdata, handles)


%= Reuse proteins from last translation spot
function button_reuse_proteins_last_Callback(hObject, eventdata, handles)

ind_trans = get(handles.listbox_translation,'Value');

if ind_trans > 1
    handles.translation(ind_trans).x_prot = handles.translation(ind_trans-1).x_prot;
    handles.translation(ind_trans).y_prot = handles.translation(ind_trans-1).y_prot;
    handles.translation(ind_trans).z_prot = handles.translation(ind_trans-1).z_prot;
    
    handles.x_prot_all = [handles.x_prot_all;handles.translation(ind_trans-1).x_prot];
    handles.y_prot_all = [handles.y_prot_all;handles.translation(ind_trans-1).y_prot];
    handles.z_prot_all = [handles.z_prot_all;handles.translation(ind_trans-1).z_prot];
    
end

 %- Save data & plot
guidata(hObject, handles);    
listbox_translation_Callback(hObject, eventdata, handles)


%= Reuse proteins from another translation spot
function button_reuse_proteins_Callback(hObject, eventdata, handles)

str_list = get(handles.listbox_translation,'String');

ind_1 = get(handles.listbox_translation,'Value');

waitfor(msgbox('Select hotspot to which proteins should be assigned and press enter.'));
ind_2 = get(handles.listbox_translation,'Value');

text_title = ['Use spots from .. ',str_list{ind_1},' .. for .. ',str_list{ind_2},'?'];

choice = questdlg(text_title, 'Reuse proteins', 'Yes','No','Yes');

if strcmp(choice,'Yes')
    
    handles.translation(ind_2).x_prot = handles.translation(ind_1).x_prot;
    handles.translation(ind_2).y_prot = handles.translation(ind_1).y_prot;
    handles.translation(ind_2).z_prot = handles.translation(ind_1).z_prot;
    
    handles.x_prot_all = [handles.x_prot_all;handles.translation(ind_1).x_prot];
    handles.y_prot_all = [handles.y_prot_all;handles.translation(ind_1).y_prot];
    handles.z_prot_all = [handles.z_prot_all;handles.translation(ind_1).z_prot];
    handles.ind_all    = [handles.ind_all;repmat(ind_2,size(handles.x_prot))];
    
    %- Save data & plot
    guidata(hObject, handles);    
    listbox_translation_Callback(hObject, eventdata, handles)
end


%= Delete last protein
function button_delete_protein_Callback(hObject, eventdata, handles)
ind_trans  = get(handles.listbox_translation,'Value');

if ~isempty(handles.translation(ind_trans).x_prot)
    handles.translation(ind_trans).x_prot(end) = [];
    handles.translation(ind_trans).y_prot(end) = [];
    handles.translation(ind_trans).z_prot(end) = [];   
    
    handles.x_prot_all(end) = [];
    handles.y_prot_all(end) = [];
    handles.z_prot_all(end) = [];
    handles.ind_all(end)    = [];
    
    
    %- Save data & plot
    guidata(hObject, handles);
    listbox_translation_Callback(hObject, eventdata, handles)
end


%= Reuse proteins from another translation spot
function button_delete_protein_all_Callback(hObject, eventdata, handles)

%- Ask user to confirm choice
choice = questdlg('Do you really want to delete ALL proteins?', mfilename, 'Yes','No','No');

if strcmp(choice,'Yes')

    ind_trans  = get(handles.listbox_translation,'Value');

    if ~isempty(handles.translation(ind_trans).x_prot)
        
        %- Delete from assosciated list
        handles.translation(ind_trans).x_prot = [];
        handles.translation(ind_trans).y_prot = [];
        handles.translation(ind_trans).z_prot = [];
               
        %- Delete from all list
        ind_delete = find(handles.ind_all == ind_trans);
        handles.x_prot_all(ind_delete) = [];
        handles.y_prot_all(ind_delete) = [];
        handles.z_prot_all(ind_delete) = [];
     
        %- Save data & plot
        guidata(hObject, handles);
        listbox_translation_Callback(hObject, eventdata, handles)
    end
end

%==========================================================================
% Analysis - hotspot quantification
%==========================================================================

%= Analyze spots
function button_analyze_Callback(hObject, eventdata, handles)

%- Loop over all translation spots
translation = handles.translation;
img_prot    = handles.img_prot;

crop = handles.img_prot.settings.avg_spots.crop;


%- Folder to save results
[dum base]  = fileparts(handles.img_prot.file_names.raw);
folder_save = fullfile(handles.img_prot.path_names.img,'_HotSpot',base);
if ~exist(folder_save); mkdir(folder_save); end

%-- Dimension of default image
dim.X = handles.img_prot.settings.avg_spots.crop.xy*2 +1;
dim.Y = handles.img_prot.settings.avg_spots.crop.xy*2 +1;
N_pix = dim.X*dim.Y;

%- Generate vectors describing the image    
[Xs,Ys] = meshgrid(dim.Y:2*dim.Y-1,dim.X:2*dim.X-1);  % Pixel-grid has on offset from 0 to allow for negative values in center position
xdata(1,:) = double(reshape(Xs,1,N_pix));
xdata(2,:) = double(reshape(Ys,1,N_pix));

options_fit.par_start.sigmax = 1; 
options_fit.par_start.sigmay = 1; 
options_fit.options = optimset('Jacobian','off','Display','off','MaxIter',100,'UseParallel','always');

bound.lb   = [0  -inf -inf  0   0  ]; 
bound.ub   = [10  inf  inf   inf inf];

%- Parameters
options_fit.pixel_size.xy = 1;
options_fit.par_start  = [];
options_fit.bound      = bound;
options_fit.fit_mode   = 'sigma_free_xz';
flag_struct.output     = 0;

par_start.sigmax = 1;

%- Integration range
x_int.min = -crop.xy; x_int.max = crop.xy;
y_int.min = -crop.xy; y_int.max = crop.xy;

for iT = 1:numel(translation)
    
   %==== Average spots
   x_prot = translation(iT).x_prot;
   y_prot = translation(iT).y_prot;
   z_prot = translation(iT).z_prot;
   
   %- Detected positions
   img_prot.cell_prop(1).spots_detected = [];
   
   img_prot.cell_prop(1).spots_detected(:,1) = y_prot;
   img_prot.cell_prop(1).spots_detected(:,2) = x_prot;
   img_prot.cell_prop(1).spots_detected(:,3) = z_prot;
   
   %- Thresholded spots
   img_prot.cell_prop(1).thresh.in      = [];
   img_prot.cell_prop(1).thresh.in      = true(size(img_prot.cell_prop(1).spots_detected(:,1)));

   [spot_avg, spot_avg_os, pixel_size_os, img_sum] = spot_3D_avg_v1(img_prot,1,[]);
   mip_avg = max(spot_avg,[],3);
   
   %==== Fit averaged spot after MIP
   
   %- Determine center of mass - starting point for center
   center_mass  = ait_centroid3d_v3(double(mip_avg),xdata);
   if not(isfield(par_start,'centerx'));  par_start.centerx = center_mass(1); end
   if not(isfield(par_start,'centery'));  par_start.centery = center_mass(2); end

   %- Min and Max of the image: starting point for amplitude and background
   img_max = max(mip_avg(:));
   img_min = (min(mip_avg(:))) * double((min(mip_avg(:))>0)) + (1*(min(mip_avg(:))<=0));

   par_start.amp = img_max-img_min; 
   par_start.bgd = img_min;  
   
   options_fit.par_start = par_start;
   
   %- Fit averaged spots
   fit_avg = spot_2D_fit_v1(mip_avg, xdata, options_fit, flag_struct);
   
   %- Calculate integrated intensity
   par_mod_int(1) = fit_avg.sigmaX; par_mod_int(2)  = fit_avg.sigmaY;   
   par_mod_int(3) = 0;              par_mod_int(4)  = 0;
   par_mod_int(5) = fit_avg.amp;    par_mod_int(6)  = 0 ;

   intint_avg = fun_Gaussian_2D_double_integral_v1(x_int,y_int,par_mod_int);
   
   %==== HOT-SPOT
   x = translation(iT).x;
   y = translation(iT).y;
   z = translation(iT).z;
   
   xmin = (x - crop.xy); xmax = (x + crop.xy);
   ymin = (y - crop.xy); ymax = (y + crop.xy);
   zmin = (z - crop.z);  zmax = (z + crop.z);
    
   %- Cropped image
   try
        hot_crop = handles.img_prot.raw(ymin:ymax,xmin:xmax,zmin:zmax);
   catch
       disp(['Hotspot - ', translation(iT).name, ' - too close to border of image. Can not be considered'])
       disp(' ')
       summary_analysis(iT,:) =    [-1  ...
                               -1  -1 -1 -1 ...
                               -1  -1 -1 -1];
         
        label_hotspots{iT,:} =      translation(iT).name; 

       continue
   end
   mip_hot  = max(hot_crop,[],3);
   
   %- Min and Max of the image: starting point for amplitude and background
   img_max   = max(mip_hot(:));
   img_min   = (min(mip_hot(:))) * double((min(mip_hot(:))>0)) + (1*(min(mip_hot(:))<=0));

   par_start.amp = img_max-img_min; 
   par_start.bgd = img_min;  
   
   options_fit.par_start = par_start;
   
   %- Fit averaged spots
   fit_hot = spot_2D_fit_v1(mip_hot,xdata,options_fit,flag_struct);
   
   %- Calculate integrated intensity
   par_mod_int(1)  = fit_hot.sigmaX;  par_mod_int(2)  = fit_hot.sigmaY;   
   par_mod_int(3)  = 0;               par_mod_int(4)  = 0;
   par_mod_int(5)  = fit_hot.amp;     par_mod_int(6)  = 0 ;

   intint_hot = fun_Gaussian_2D_double_integral_v1(x_int,y_int,par_mod_int);
   
   %- Summary
   summary_analysis(iT,:) =    [x y z... 
                               intint_hot/intint_avg  ...
                               intint_avg  fit_avg.sigmaX fit_avg.amp fit_avg.bgd ...
                               intint_hot  fit_hot.sigmaX fit_hot.amp fit_hot.bgd];
   
                           
   label_hotspots{iT,:} = translation(iT).name; 
                           
   %- Plot results
   figure, set(gcf, 'visible','off')
   subplot(1,2,1)
   imshow(mip_hot,[])
   colorbar
   title(handles.img_prot.file_names.raw, 'Interpreter', 'none')
    
   subplot(1,2,2)
   imshow(mip_avg,[])
   colorbar
   title('Averaged protein')
   
   name_save = fullfile(folder_save,translation(iT).name);
   saveas(gcf,name_save,'png')
   close(gcf)

end

%- Save summary
handles.summary_analysis = summary_analysis;
handles.label_hotspots   = label_hotspots;
handles.folder_save = folder_save;
guidata(hObject, handles);


%==========================================================================
% Analysis - CoLoc
%==========================================================================


%- New cell
function button_cl_cell_new_Callback(hObject, eventdata, handles)

param.reg_type = 'Freehand';
param.h_axes   = handles.axes_protein;
param.pos      = [];

reg_result = FQ_draw_region_v1(param);

position = reg_result.position;

if ~isempty(position)         

    %- Get current list
    str_list = get(handles.listbox_cl_cells,'String');
    N_Cell   = size(str_list,1);
    ind_cell = N_Cell+1;

    %- Save position
    handles.cell_prop(ind_cell).x = round(position(:,1))';  % v3: Has to be a row vector to agree with read-in from files
    handles.cell_prop(ind_cell).y = round(position(:,2))';  % v3: Has to be a row vector to agree with read-in from files

    %- Add entry at the end and update list
    str_cell = ['Cell_', num2str(handles.cell_counter )];
    str_list{ind_cell} = str_cell;

    set(handles.listbox_cl_cells,'String',str_list)
    set(handles.listbox_cl_cells,'Value',ind_cell)

    handles.cell_prop(ind_cell).name = str_cell;
    handles.cell_counter = handles.cell_counter+1;
    
    %- Show plot
    listbox_cl_cells_Callback(hObject, eventdata, handles)
    
    %- Save results
    guidata(hObject, handles);
end


%- Delete cell
function button_cl_cell_delete_Callback(hObject, eventdata, handles)
%- Ask user to confirm choice
choice = questdlg('Do you really want to delete this cell?', 'HotSpot', 'Yes','No','No');

if strcmp(choice,'Yes')
    
    %- Extract index of highlighted cell
    str_list = get(handles.listbox_cl_cells,'String');
    ind_sel  = get(handles.listbox_cl_cells,'Value');
    
    %- Delete highlighted cell
    str_list(ind_sel) = [];
    
    handles.cell_prop(ind_sel) = [];   
    set(handles.listbox_cl_cells,'String',str_list)
    
    %- Make sure that pointer is not outside of defined cells
    N_str = length(str_list);
    if ind_sel > N_str     
        set(handles.listbox_cl_cells,'Value',N_str)
    end
    
    %- Show plot
    listbox_cl_cells_Callback(hObject, eventdata, handles) 
    
    %- Save results
    guidata(hObject, handles);
end


%- Listbox with cells
function listbox_cl_cells_Callback(hObject, eventdata, handles)
handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);


%- Add proteins
function button_cl_protein_add_Callback(hObject, eventdata, handles)

%- Check if cell is defined
ind_cell = get(handles.listbox_cl_cells,'Value');
if isempty(handles.cell_prop(ind_cell).x)
    errordlg('No cell defined!','HotSpot')
    return
end

%- Define proteins
x = [];
y = [];

while 1
    [x_loop,y_loop,button,axn] = ginputax(handles.axes_protein,1);
   
    %- Empty position - button pressed
    if isempty(x_loop); break; end
    
    %- Check if in actual axes (and not in axes for RNA)
    if isempty(axn); break; end
    
    %- Make sure not outside of the axes but in the GUI
    if  (x_loop < 0.5  ||  x_loop > (handles.img_prot.dim.X+0.5) && ...
       y_loop < 0.5  ||  y_loop > (handles.img_prot.dim.Y+0.5))
        
        break
    end
    
    %- Save positions
    x=[x;x_loop]; y=[y;y_loop];

    %- Show image
    hold on
        plot(handles.axes_protein,x_loop,y_loop,'g+')
    hold off
 
end

%- Return if no position was defined
if isempty(x); return; end

%- Analyze positions
x = round(x);
y = round(y); 

%- Go over spots and find maximum
crop = handles.img_prot.settings.avg_spots.crop;

for iS = 1:numel(x)
    
    %- Cropping window
    xmin = (x(iS) - crop.xy); 
    xmax = (x(iS) + crop.xy);
    
    ymin = (y(iS) - crop.xy);
    ymax = (y(iS) + crop.xy);
    
    %- Cropped image
    try
        img_crop = handles.img_prot.raw_proj_z(ymin:ymax,xmin:xmax,:);
    catch
        listbox_translation_Callback(hObject, eventdata, handles)
        errordlg('Spot close to edge of image (in X,Y, or Z)','HotSpot')
        return
    end
    
    %figure, imshow(max(img_crop,[],3),[])
    
    %- Find maximum position in this image
    [dum, Ixyz] = max(img_crop(:));
    [Iy,Ix,Iz]  = ind2sub(size(img_crop),Ixyz);
    
    %- Correct original estimates by this offset
    x_new(iS,1) = x(iS) + (Ix-crop.xy-1);
    y_new(iS,1) = y(iS) + (Iy-crop.xy-1);
    
end


if isempty(handles.cell_prop(ind_cell).x_prot)
    handles.cell_prop(ind_cell).x_prot = round(x_new);
    handles.cell_prop(ind_cell).y_prot = round(y_new);
else
    handles.cell_prop(ind_cell).x_prot = [handles.cell_prop(ind_cell).x_prot;round(x_new)];
    handles.cell_prop(ind_cell).y_prot = [handles.cell_prop(ind_cell).y_prot;round(y_new)];
end

%- Save data & plot
guidata(hObject, handles);
listbox_cl_cells_Callback(hObject, eventdata, handles)


%=== ADD rna
function button_cl_rna_add_Callback(hObject, eventdata, handles)

%- Check if cell is defined
ind_cell = get(handles.listbox_cl_cells,'Value');
if isempty(handles.cell_prop(ind_cell).x)
    errordlg('No cell defined!','HotSpot')
    return
end

%- Define proteins
x = [];
y = [];

while 1
    [x_loop,y_loop,button,axn] = ginputax(handles.axes_rna,1);
   
    %- Empty position - button pressed
    if isempty(x_loop); break; end
    
    %- Check if in actual axes (and not in axes for RNA)
    if isempty(axn); break; end
    
    %- Make sure not outside of the axes but in the GUI
    if  (x_loop < 0.5  ||  x_loop > (handles.img_prot.dim.X+0.5) && ...
       y_loop < 0.5  ||  y_loop > (handles.img_prot.dim.Y+0.5))
        
        break
    end
    
    %- Save positions
    x=[x;x_loop]; y=[y;y_loop];

    %- Show image
    hold on
        plot(handles.axes_rna,x_loop,y_loop,'xr')
    hold off
 
end

%- Return if no position was defined
if isempty(x); return; end

%- Analyze positions
x = round(x);
y = round(y); 

%- Go over spots and find maximum
crop = handles.img_prot.settings.avg_spots.crop;

for iS = 1:numel(x)
    
    %- Cropping window
    xmin = (x(iS) - crop.xy); 
    xmax = (x(iS) + crop.xy);
    
    ymin = (y(iS) - crop.xy);
    ymax = (y(iS) + crop.xy);
    
    %- Cropped image
    try
        img_crop = handles.img_rna.raw_proj_z(ymin:ymax,xmin:xmax,:);
    catch
        listbox_translation_Callback(hObject, eventdata, handles)
        errordlg('Spot close to edge of image (in X,Y, or Z)','HotSpot')
        return
    end
    
    %figure, imshow(max(img_crop,[],3),[])
    
    %- Find maximum position in this image
    [dum, Ixyz] = max(img_crop(:));
    [Iy,Ix,Iz]  = ind2sub(size(img_crop),Ixyz);
    
    %- Correct original estimates by this offset
    x_new(iS,1) = x(iS) + (Ix-crop.xy-1);
    y_new(iS,1) = y(iS) + (Iy-crop.xy-1);
    
end


if isempty(handles.cell_prop(ind_cell).x_rna)
    handles.cell_prop(ind_cell).x_rna = round(x_new);
    handles.cell_prop(ind_cell).y_rna = round(y_new);
else
    handles.cell_prop(ind_cell).x_rna = [handles.cell_prop(ind_cell).x_rna;round(x_new)];
    handles.cell_prop(ind_cell).y_rna = [handles.cell_prop(ind_cell).y_rna;round(y_new)];
end

%- Save data & plot
guidata(hObject, handles);
listbox_cl_cells_Callback(hObject, eventdata, handles)


%- Delete last protein
function button_cl_protein_delete_Callback(hObject, eventdata, handles)

ind_sel  = get(handles.listbox_cl_cells,'Value');

if ~isempty(handles.cell_prop(ind_sel).x_prot)
    handles.cell_prop(ind_sel).x_prot(end) = [];
    handles.cell_prop(ind_sel).y_prot(end) = [];
end

%- Save data & plot
guidata(hObject, handles);
listbox_cl_cells_Callback(hObject, eventdata, handles)

    
%- Delete all proteins from current cell
function button_cl_protein_delete_all_Callback(hObject, eventdata, handles)
 
choice = questdlg('Do you really want to delete ALL proteins?', mfilename, 'Yes','No','No');

if strcmp(choice,'Yes')
    
    ind_sel  = get(handles.listbox_cl_cells,'Value');
    if ~isempty(handles.cell_prop(ind_sel).x_prot)
        handles.cell_prop(ind_sel).x_prot = [];
        handles.cell_prop(ind_sel).y_prot = [];
    end

    %- Save data & plot
    guidata(hObject, handles);
    listbox_cl_cells_Callback(hObject, eventdata, handles)
end
    

%- Delete last RNA
function button_cl_rna_delete_Callback(hObject, eventdata, handles)

ind_sel  = get(handles.listbox_cl_cells,'Value');

if ~isempty(handles.cell_prop(ind_sel).x_rna)
    handles.cell_prop(ind_sel).x_rna(end) = [];
    handles.cell_prop(ind_sel).y_rna(end) = [];
end

%- Save data & plot
guidata(hObject, handles);
listbox_cl_cells_Callback(hObject, eventdata, handles)


%- Delete all RNA
function button_cl_rna_delete_all_Callback(hObject, eventdata, handles)
choice = questdlg('Do you really want to delete ALL RNAs?', mfilename, 'Yes','No','No');

if strcmp(choice,'Yes')
    ind_sel  = get(handles.listbox_cl_cells,'Value');
    if ~isempty(handles.cell_prop(ind_sel).x_rna)
        handles.cell_prop(ind_sel).x_rna = [];
        handles.cell_prop(ind_sel).y_rna = [];
    end

    %- Save data & plot
    guidata(hObject, handles);
    listbox_cl_cells_Callback(hObject, eventdata, handles)
end
    



%% ==========================================================================
%   Saving and loading - hotspots
%  ==========================================================================

%=== SAVE SPOTS
function menu_save_spots_Callback(hObject, eventdata, handles)

%- File-name
[dum base] = fileparts(handles.img_prot.file_names.raw);
file_name_default = fullfile(handles.img_prot.path_names.img,[base,'__HOTspots_SPOTS.txt']);

[file_save,path_save] = uiputfile(file_name_default,'Save defined spots of HOTspots analysis');
if file_save == 0; return; end  
file_name_full = fullfile(path_save,file_save);

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'HotSpots analysis - SPOTS\t %s \n\n', date);
fprintf(fid,'ImgPROT\t %s \n', handles.img_prot.file_names.raw);
fprintf(fid,'ImgRNA\t %s \n', handles.img_rna.file_names.raw);
fprintf(fid,'PATH\t %s \n', handles.img_prot.path_names.img);

translation= handles.translation;

for iT = 1:numel(translation)

   %==== Location of hot-spot
   fprintf(fid,'HotSpotName\t%s \n X\tY\tZ \n', translation(iT).name);
   fprintf(fid,'%g\t%g\t%g\t\n',translation(iT).x, translation(iT).y,translation(iT).z);
         
   %==== Location & intensity of proteins
   idx = sub2ind(size(handles.img_prot.raw), ...
                      translation(iT).y_prot, ...
                      translation(iT).x_prot, ...
                      translation(iT).z_prot);
   prot_int = handles.img_prot.raw(idx);
   
   
   fprintf(fid,'Proteins\n');
   N_par = numel( translation(iT).y_prot);
   string_write = ['%s',repmat('\t%g',1,N_par),'\n'];
   
   fprintf(fid,string_write, 'X',translation(iT).x_prot);
   fprintf(fid,string_write, 'Y',translation(iT).y_prot);
   fprintf(fid,string_write, 'Z',translation(iT).z_prot);
  fprintf(fid,string_write, 'INT', prot_int);
    
end
    

%=== Load SPOTS
function menu_load_spots_Callback(hObject, eventdata, handles)

%- Load file
[file_name,path_name] = uigetfile({'*.txt'},'HotSpot - select file with spot localizations');    
if file_name == 0; return; end

%- Open file
fid  = fopen(fullfile(path_name,file_name),'r');

if fid == -1
    warndlg('File cannot be opened','HotSpot'); 
    disp(['File : ' , file_name])
    disp(['Path : ' , path_name])
    return
else
    
    disp('=== HotSpot - OPENING SPOT FILE')
    disp(['File : ' , file_name])
    disp(['Path : ' , path_name])
    
    handles.translation = struct('name','','x','','y','','z','','x_prot','','y_prot','','z_prot','');
    ind_translation     = 0;
    
    %- Go over lines
    tline = fgetl(fid);
    
	while ischar(tline)
    
        disp(tline)
        
        if not(isempty(strfind(tline, 'HotSpotName')))
                
            %=== Change index
            ind_translation = ind_translation + 1;
            
            %=== hotspot: name
            k = strfind(tline, sprintf('\t') );
            handles.translation(ind_translation).name = tline(k+1:end);  
            
            str_list{ind_translation} = handles.translation(ind_translation).name;
            
            %-  hotspot: XYZ-coordinates
            fgetl(fid);     % Identifier X,Y,Z
            tline = fgetl(fid);
            pos   = str2num(tline);  
            handles.translation(ind_translation).x = pos(1);
            handles.translation(ind_translation).y = pos(2);
            handles.translation(ind_translation).z = pos(3);
            
            %- Identifier protein
            fgetl(fid);
            tline = fgetl(fid);
            k     = [strfind(tline, sprintf('\t') ),length(tline)];
            pos = str2num(tline(k(1)+1:k(end)));
            handles.translation(ind_translation).x_prot = pos';
            
            tline = fgetl(fid);
            k     = [strfind(tline, sprintf('\t') ),length(tline)];
            pos = str2num(tline(k(1)+1:k(end)));
            handles.translation(ind_translation).y_prot = pos';
            
            tline = fgetl(fid);
            k     = [strfind(tline, sprintf('\t') ),length(tline)];
            pos = str2num(tline(k(1)+1:k(end)));
            handles.translation(ind_translation).z_prot = pos';           

        end
        
        %- Get next line
        tline = fgetl(fid);
        
        % Read one more line if this is a newer file containing also the intensity values
        if ischar(tline) && contains(tline,'INT')
            tline = fgetl(fid);
        end
    end    
end
   
if isempty(str_list)
    warndlg('No HotSpots found. Are you sure this is a valid file?','Load saved HotSpots')
else 
    %-- save and populate GUI
    guidata(hObject, handles);
    set(handles.listbox_translation ,'String',str_list)
    set(handles.listbox_translation ,'Value',1)

    listbox_translation_Callback(hObject, eventdata, handles) 
end


%=== SAVE RESULTS
function menu_save_results_Callback(hObject, eventdata, handles)

%- File-name
[dum base] = fileparts(handles.img_prot.file_names.raw);
file_name_default = fullfile(handles.img_prot.path_names.img,[base,'__HOTspots_RESULTS.txt']);

[file_save,path_save] = uiputfile(file_name_default,'Save results of HOTspots analysis');
if file_save == 0; return; end  
file_name_full = fullfile(path_save,file_save);

%- Prepare data
N_par        = size(handles.summary_analysis,2);
cell_data    = num2cell(handles.summary_analysis);  
cell_write   = [handles.label_hotspots,cell_data];      
cell_write   = cell_write';   %- fprintf works on colums - data has to therefore be transformed 
string_write = ['%s',repmat('\t%g',1,N_par),'\n'];

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'HotSpots analysis\t %s \n', date);
fprintf(fid,'Name\tx_pos\ty_pos\tz_pos\tNprot\tAvg_Intint\tAvg_sigma\tAvg_amp\tAvg_bgd\tHot_Intint\tHot_sigma\tHot_amp\tHot_bgd\n');
fprintf(fid,string_write,cell_write{:});             
fclose(fid);     


%% ==========================================================================
%   Saving and loading - CoLoc
%  ==========================================================================

%- Save summary
function menu_cl_save_summary_Callback(hObject, eventdata, handles)

%- File-name
img_name = handles.img_prot.file_names.raw;
[dum base] = fileparts(handles.img_prot.file_names.raw);
file_name_default = fullfile(handles.img_prot.path_names.img,[base,'__CoLoc_SUMMARY.txt']);

[file_save,path_save] = uiputfile(file_name_default,'Save results of CoLoc analysis');
if file_save == 0; return; end  
file_name_full = fullfile(path_save,file_save);

cell_prop = handles.cell_prop;

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'CoLoc analysis - Summary\t %s \n', date);
fprintf(fid,'Image-Name\tCell-Name\tN_Col\tN_NotCol\n');
for iC = 1:numel(cell_prop)
    N_Col    = numel(cell_prop(iC).x_prot);
    N_NotCol = numel(cell_prop(iC).x_rna);
    fprintf(fid,'%s\t%s\t%g\t%g\n',img_name,cell_prop(iC).name,N_Col,N_NotCol);    
end
fclose(fid);     


%- Save spots
function menu_cl_save_spots_Callback(hObject, eventdata, handles)
%- File-name
img_name = handles.img_prot.file_names.raw;
[dum base] = fileparts(handles.img_prot.file_names.raw);
file_name_default = fullfile(handles.img_prot.path_names.img,[base,'__CoLoc_SPOTS.txt']);

[file_save,path_save] = uiputfile(file_name_default,'Save results of CoLoc analysis');
if file_save == 0; return; end  
file_name_full = fullfile(path_save,file_save);

cell_prop = handles.cell_prop;

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'CoLoc analysis - Summary\t %s \n', date);
fprintf(fid,'Image-Name\t%s\n',img_name);
for iC = 1:numel(cell_prop)
  
    fprintf(fid,'Cell-Name\t%s\n',cell_prop(iC).name);  
    
    %- Cell position
    N_par = numel(cell_prop(iC).x);
    string_write = ['%s\t%g',repmat('\t%g',1,N_par-1),'\n'];
    fprintf(fid,string_write, 'CELL_X',cell_prop(iC).x);
    fprintf(fid,string_write, 'CELL_Y',cell_prop(iC).y);
    
    %- Protein position
    N_par = numel(cell_prop(iC).x_prot);
    string_write = ['%s\t%g',repmat('\t%g',1,N_par-1),'\n'];
    fprintf(fid,string_write, 'PROT_X',cell_prop(iC).x_prot);
    fprintf(fid,string_write, 'PROT_Y',cell_prop(iC).y_prot);
    
    %- RNA position
    N_par = numel(cell_prop(iC).x_rna);
    string_write = ['%s\t%g',repmat('\t%g',1,N_par-1),'\n'];
    fprintf(fid,string_write, 'RNA_X',cell_prop(iC).x_rna);
    fprintf(fid,string_write, 'RNA_Y',cell_prop(iC).y_rna);
    
end
fclose(fid);     


%==========================================================================
% Plotting
%==========================================================================\

%== Plot image -  hotspot analysis
function handles = plot_image(hObject, eventdata, handles)
if get(handles.popup_analysis,'Value') == 1
    handles = plot_image_hs(hObject, eventdata, handles);
elseif get(handles.popup_analysis,'Value') == 2
    handles = plot_image_cl(hObject, eventdata, handles);
end


%== Plot results from coloc analysis
function handles = plot_image_cl(hObject, eventdata, handles)

axes(handles.axes_protein); 
v = axis(handles.axes_protein);


ind_sel  = get(handles.listbox_cl_cells,'Value');
if isempty(handles.cell_prop)
    x_cell = []; y_cell = [];
    x_prot = []; y_prot = [];
    x_rna = [];  y_rna = [];
else
    x_cell = handles.cell_prop(ind_sel).x;
    y_cell = handles.cell_prop(ind_sel).y;
    
    x_prot = handles.cell_prop(ind_sel).x_prot;
    y_prot = handles.cell_prop(ind_sel).y_prot;

    x_rna = handles.cell_prop(ind_sel).x_rna;
    y_rna = handles.cell_prop(ind_sel).y_rna;
    
end


%== Plot Protein data
if ~isempty(handles.img_prot.raw)
    
    %- Contrast
    slider_prot_min = get(handles.slider_contrast_prot_min,'Value');
    slider_prot_max = get(handles.slider_contrast_prot_max,'Value');

    Im_prot_min = slider_prot_min*handles.img_prot_diff+handles.img_prot_min;
    Im_prot_max = slider_prot_max*handles.img_prot_diff+handles.img_prot_min;

    if Im_prot_max < Im_prot_min; Im_prot_max = Im_prot_min+1; end
    
    img_plot =  handles.img_prot.raw_proj_z;

    %- Show image
    axes(handles.axes_protein); 
    imshow(img_plot,[Im_prot_min Im_prot_max]); 
    
    hold on
        for iC=1:numel(handles.cell_prop)
            x_dum = handles.cell_prop(iC).x;
            y_dum = handles.cell_prop(iC).y;
            plot(x_dum,y_dum,'-y'); 
        end
    
        plot(x_cell,y_cell,'-b');    
        plot(x_prot,y_prot,'og');
        plot(x_rna,y_rna,'or');
    hold off
end

%== Plot RNA data
if ~isempty(handles.img_rna.raw)
    
    %- Contrast
    slider_rna_min = get(handles.slider_contrast_rna_min,'Value');
    slider_rna_max = get(handles.slider_contrast_rna_max,'Value');

    Im_rna_min = slider_rna_min*handles.img_rna_diff+handles.img_rna_min;
    Im_rna_max = slider_rna_max*handles.img_rna_diff+handles.img_rna_min;

    if Im_rna_max < Im_rna_min; Im_rna_max = Im_rna_min+1; end
    
     img_plot =  handles.img_rna.raw_proj_z;

    %- Show image
    axes(handles.axes_rna); 
    cla   
    imshow(img_plot,[Im_rna_min Im_rna_max])   

    hold on
       for iC=1:numel(handles.cell_prop)
            x_dum = handles.cell_prop(iC).x;
            y_dum = handles.cell_prop(iC).y;
            plot(x_dum,y_dum,'-y'); 
        end
    
        plot(x_cell,y_cell,'b');
        plot(x_prot,y_prot,'og');
        plot(x_rna,y_rna,'or');
    hold off
    
    %- Same zoom as before
    if not(handles.status_first_plot)
        axis(v);
    end

end

%== Plot results from hotspot analysis
function handles = plot_image_hs(hObject, eventdata, handles)

%- Get track and current position
ind_frame = str2double(get(handles.text_frame,'String'));
show_proteins = get(handles.checkbox_show_proteins,'Value');
axes(handles.axes_protein); 
v = axis(handles.axes_protein);

%== Plot translational spot and associated proteins
x_trans=[];       y_trans=[];
x_trans_all = []; y_trans_all = [];
x_prot = [];      y_prot = [];
    
ind_sel  = get(handles.listbox_translation,'Value');
% checkbox_show_proteins
if ~isempty(handles.translation)
    if ~isempty(handles.translation(ind_sel).x)

        x_trans = handles.translation(ind_sel).x;
        y_trans = handles.translation(ind_sel).y;

        x_trans_all = [handles.translation.x];
        y_trans_all = [handles.translation.y];
        
        
        if show_proteins && ~isempty(handles.translation(ind_sel).x_prot)
            x_prot = handles.translation(ind_sel).x_prot;
            y_prot = handles.translation(ind_sel).y_prot;
        end

    end
end

%== Plot Protein data
if ~isempty(handles.img_prot.raw)
    
    %- Contrast
    slider_prot_min = get(handles.slider_contrast_prot_min,'Value');
    slider_prot_max = get(handles.slider_contrast_prot_max,'Value');

    Im_prot_min = slider_prot_min*handles.img_prot_diff+handles.img_prot_min;
    Im_prot_max = slider_prot_max*handles.img_prot_diff+handles.img_prot_min;

    if Im_prot_max < Im_prot_min; Im_prot_max = Im_prot_min+1; end
    
    img_plot =  handles.img_prot.raw(:,:,ind_frame);

    %- Show image
    axes(handles.axes_protein); 
    imshow(img_plot,[Im_prot_min Im_prot_max]);   

    hold on
        
        plot(x_trans_all,y_trans_all,'oy')      
        plot(x_trans,y_trans,'or')
        
        %plot(handles.x_prot_all,handles.y_prot_all,'y+')
        if show_proteins && ~isempty(handles.translation(ind_sel).x_prot)
            plot(x_prot,y_prot,'g+')
        end
        
    hold off
    
    %- Same zoom as before
    if not(handles.status_first_plot)
        axis(v);
    end
end


%== Plot RNA data
if ~isempty(handles.img_rna.raw)
    
    %- Contrast
    slider_rna_min = get(handles.slider_contrast_rna_min,'Value');
    slider_rna_max = get(handles.slider_contrast_rna_max,'Value');

    Im_rna_min = slider_rna_min*handles.img_rna_diff+handles.img_rna_min;
    Im_rna_max = slider_rna_max*handles.img_rna_diff+handles.img_rna_min;

    if Im_rna_max < Im_rna_min; Im_rna_max = Im_rna_min+1; end
    
    img_plot =  handles.img_rna.raw(:,:,ind_frame);

    %- Show image
    axes(handles.axes_rna); 
    cla   
    imshow(img_plot,[Im_rna_min Im_rna_max])   
    
    hold on
        plot(x_trans_all,y_trans_all,'oy')      
        plot(x_trans,y_trans,'or')
         
        % plot(handles.x_prot_all,handles.y_prot_all,'y+')
         plot(x_prot,y_prot,'g+')
         
         
    hold off
    
    %- Same zoom as before
    if not(handles.status_first_plot)
        axis(v);
    end

end

%- Save everything
handles.status_first_plot = 0;

%- Save data
guidata(hObject, handles);


%== Zoom
function button_zoom_Callback(hObject, eventdata, handles)
if handles.status_zoom == 0
    h_zoom = zoom;
    set(h_zoom,'Enable','on');
    handles.status_zoom = 1;
    handles.status_pan  = 0;
    handles.h_zoom      = h_zoom;
else
    set(handles.h_zoom,'Enable','off');    
    handles.status_zoom = 0;
end
guidata(hObject, handles);


%== Pan
function button_pan_Callback(hObject, eventdata, handles)
if handles.status_pan == 0
    h_pan = pan;
    set(h_pan,'Enable','on');
    handles.status_pan  = 1;
    handles.status_zoom = 0;
    handles.h_pan      = h_pan;    
else
    set(handles.h_pan,'Enable','off');    
    handles.status_pan = 0;
end
guidata(hObject, handles);


%=== Listbox with translational spots
function listbox_translation_Callback(hObject, eventdata, handles)
handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);


%% ==========================================================================
%   Control for frames
%  ==========================================================================

%== Slider for slice
function slider_slice_Callback(hObject, eventdata, handles)
slider_value = get(handles.slider_slice,'Value');
ind_slice = 1+round(slider_value*(handles.img_prot_NZ-1));
set(handles.text_frame,'String',num2str(ind_slice));

handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

%== Up one slice
function button_slice_incr_Callback(hObject, eventdata, handles)

%- Check how many slices
N_slice      = handles.img_prot_NZ;
track_start  = 1;

%- Get next value for slice
ind_slice = str2double(get(handles.text_frame,'String'))+1;
if ind_slice > N_slice;ind_slice = N_slice;end
set(handles.text_frame,'String',ind_slice);

%-Update slider
slider_value = double(ind_slice-track_start)/double((N_slice-1));
set(handles.slider_slice,'Value',slider_value);

%- Save and plot image
handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

%== Down one slice
function button_slice_decr_Callback(hObject, eventdata, handles)

%- Get next value for slice
ind_slice = str2double(get(handles.text_frame,'String'))-1;
if ind_slice <1;ind_slice = 1;end
set(handles.text_frame,'String',ind_slice);

%-Update slider
slider_value = double(ind_slice-1)/double(handles.img_prot_NZ-1);
set(handles.slider_slice,'Value',slider_value);

%- Save and plot image
handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

  

%% ==========================================================================
%   Control for contrast
%  ==========================================================================

%== Slider constrast min: rna
function slider_contrast_rna_min_Callback(hObject, eventdata, handles)
slider_min = get(handles.slider_contrast_rna_min,'Value');

img_min  = handles.img_rna_min;
img_diff = handles.img_rna_diff;

contr_min = slider_min*img_diff+img_min;
set(handles.text_contr_rna_min,'String',num2str(round(contr_min)));

handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

%== Slider constrast max: rna
function slider_contrast_rna_max_Callback(hObject, eventdata, handles)
slider_min = get(handles.slider_contrast_rna_min,'Value');
slider_max = get(handles.slider_contrast_rna_max,'Value');

img_min  = handles.img_rna_min;
img_diff = handles.img_rna_diff;

contr_min = slider_min*img_diff+img_min;
contr_max = slider_max*img_diff+img_min;

if contr_max < contr_min
    contr_max = contr_min+1;
end
set(handles.text_contr_rna_max,'String',num2str(round(contr_max)));

handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

%== Slider constrast min: protein
function slider_contrast_prot_min_Callback(hObject, eventdata, handles)
slider_min = get(handles.slider_contrast_prot_min,'Value');

img_min  = handles.img_prot_min;
img_diff = handles.img_prot_diff;

contr_min = slider_min*img_diff+img_min;
set(handles.text_contr_prot_min,'String',num2str(round(contr_min)));

handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

%== Slider constrast max: protein
function slider_contrast_prot_max_Callback(hObject, eventdata, handles)
slider_min = get(handles.slider_contrast_prot_min,'Value');
slider_max = get(handles.slider_contrast_prot_max,'Value');

img_min  = handles.img_prot_min;
img_diff = handles.img_prot_diff;

contr_min = slider_min*img_diff+img_min;
contr_max = slider_max*img_diff+img_min;

if contr_max < contr_min
    contr_max = contr_min+1;
end
set(handles.text_contr_prot_max,'String',num2str(round(contr_max)));

handles = plot_image(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in checkbox_show_proteins.
function checkbox_show_proteins_Callback(hObject, eventdata, handles)
handles = plot_image_hs(hObject, eventdata, handles);
guidata(hObject, handles);

%==========================================================================
% Not used
%==========================================================================\

function slider_contrast_rna_max_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_contrast_rna_min_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_contrast_prot_max_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_contrast_prot_min_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_slice_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function text_frame_Callback(hObject, eventdata, handles)

function text_frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_translation_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_analysis_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function menu_coloc_Callback(hObject, eventdata, handles)

function listbox_cl_cells_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



