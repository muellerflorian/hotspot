%% Specify file
img_prot    =  FQ_img;
status_file = img_prot.load_img([],'raw');

if ~status_file; return; end
cd(img_prot.path_names.img)    

%% Load coordinates
[file_coord, path_coord] = uigetfile('.*csv','Specify file with coordinates',mfilename);

if file_coord == 0; return; end

%- Open file
fid  = fopen(fullfile(path_coord,file_coord),'r');

if fid == -1
    warndlg('File cannot be opened',mfilename); 
    disp(['File   : ' , file_coord])
    disp(['Folder : ' , path_coord])
    fclose(fid);
    return
else
    
    %- Read in as text
    C = textscan(fid,'%s%s%s%s','Delimiter',';','HeaderLines',1,'CollectOutput',true);
    
    %- Get data, replace NA by NaN
    data_det = C{1};
    data_det = strrep(data_det,'NA','NaN');
    
    data_det = str2double(data_det);
    fclose(fid);
end




%% Preparation for fit

img_prot.settings.avg_spots.crop.xy = 3;

%-- Dimension of default image
dim.X = img_prot.settings.avg_spots.crop.xy*2 +1;
dim.Y = img_prot.settings.avg_spots.crop.xy*2 +1;
N_pix = dim.X*dim.Y;

%- Generate vectors describing the image    
[Xs,Ys] = meshgrid(dim.Y:2*dim.Y-1,dim.X:2*dim.X-1);  % Pixel-grid has on offset from 0 to allow for negative values in center position
xdata_full(1,:) = double(reshape(Xs,1,N_pix));
xdata_full(2,:) = double(reshape(Ys,1,N_pix));

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
crop = img_prot.settings.avg_spots.crop;
x_int.min = -crop.xy; x_int.max = crop.xy;
y_int.min = -crop.xy; y_int.max = crop.xy;



%% Loop over all detections - free fit
 
summary_quant = [];

disp('FITTING DATA - free fitting parameters....')

for iD=1:size(data_det,1)
    
    if any(isnan(data_det(iD,:)))
        summary_quant(iD,1:9) = -1;
    else
    
        ind = data_det(iD,4);

        x   = data_det(iD,1)+1;
        y   = data_det(iD,2)+1;        

       xmin = (x - crop.xy); xmax = (x + crop.xy);
       ymin = (y - crop.xy); ymax = (y + crop.xy);

       %- Cropped image
       status_changed = 1;
       if ymin<1;ymin=1; end
       if xmin<1;xmin=1; status_changed=1;end
       if xmax>img_prot.dim.X;xmax=img_prot.dim.X; status_changed=1;end
       if ymax>img_prot.dim.Y;xmay=img_prot.dim.Y; status_changed=1;end

       if ~status_changed
           xdata = xdata_full;
       else
           xdata = [];
       end

       %- Cropped image
       img_crop = img_prot.raw(ymin:ymax,xmin:xmax,ind);
       %figure, imshow(img_crop,[])

       %- Min and Max of the image: starting point for amplitude and background
       img_max   = max(img_crop(:));
       img_min   = (min(img_crop(:))) * double((min(img_crop(:))>0)) + (1*(min(img_crop(:))<=0));

       par_start.amp = img_max-img_min; 
       par_start.bgd = img_min;  

       options_fit.par_start = par_start;

       %- Fit averaged spots
       flag_struct.output = 0;
       fit_spot = spot_2D_fit_v1(img_crop,xdata,options_fit,flag_struct);

       %- Calculate integrated intensity
       par_mod_int(1)  = fit_spot.sigmaX;  par_mod_int(2)  = fit_spot.sigmaY;   
       par_mod_int(3)  = 0;               par_mod_int(4)  = 0;
       par_mod_int(5)  = fit_spot.amp;     par_mod_int(6)  = 0 ;

       intint_hot = fun_Gaussian_2D_double_integral_v1(x_int,y_int,par_mod_int);

       %- Summarize
       summary_quant(iD,1:4) = data_det(iD,:);
       summary_quant(iD,5)    = intint_hot;
       summary_quant(iD,6)    = fit_spot.sigmaX;
       summary_quant(iD,7)    = fit_spot.amp;
       summary_quant(iD,8)    = fit_spot.bgd;
       summary_quant(iD,9)    = max(img_crop(:));
    end
end
disp(' .... DONE!')


%% Restrict analysis
fit_limits.sigma_xy_min = median(summary_quant(:,6)) - 0.5*std(summary_quant(:,6));
fit_limits.sigma_xy_max = median(summary_quant(:,6)) + 0.5*std(summary_quant(:,6));
fit_limits.bgd_min      = median(summary_quant(:,8));
fit_limits.bgd_max      = median(summary_quant(:,8))     + 1;    

options_fit.bound.lb   = [fit_limits.sigma_xy_min  -inf -inf  0   fit_limits.bgd_min  ]; 
options_fit.bound.ub   = [fit_limits.sigma_xy_max inf  inf   inf  fit_limits.bgd_max];

%== Loop over all
disp('FITTING DATA - RESTRICTED FITTING PARAMETERS ....')
summary_quant_restrict = [];
for iD=1:size(data_det,1)
    
     if any(isnan(data_det(iD,:)))
        summary_quant_restrict(iD,1:9) = -1;
    else
    
        ind = data_det(iD,4);

        x   = data_det(iD,1)+1;
        y   = data_det(iD,2)+1;        

       xmin = (x - crop.xy); xmax = (x + crop.xy);
       ymin = (y - crop.xy); ymax = (y + crop.xy);

       %- Cropped image
       status_changed = 1;
       if ymin<1;ymin=1; end
       if xmin<1;xmin=1; status_changed=1;end
       if xmax>img_prot.dim.X;xmax=img_prot.dim.X; status_changed=1;end
       if ymax>img_prot.dim.Y;xmay=img_prot.dim.Y; status_changed=1;end

       if ~status_changed
           xdata = xdata_full;
       else
           xdata = [];
       end

       img_crop = img_prot.raw(ymin:ymax,xmin:xmax,ind);
       %figure, imshow(img_crop,[])

       %- Min and Max of the image: starting point for amplitude and background
       img_max   = max(img_crop(:));
       img_min   = (min(img_crop(:))) * double((min(img_crop(:))>0)) + (1*(min(img_crop(:))<=0));

       par_start.amp = img_max-img_min; 
       par_start.bgd = img_min;  

       options_fit.par_start = par_start;

       %- Fit averaged spots
       flag_struct.output = 0;
       fit_spot = spot_2D_fit_v1(img_crop,xdata,options_fit,flag_struct);

       %- Calculate integrated intensity
       par_mod_int(1)  = fit_spot.sigmaX;  par_mod_int(2)  = fit_spot.sigmaY;   
       par_mod_int(3)  = 0;               par_mod_int(4)  = 0;
       par_mod_int(5)  = fit_spot.amp;     par_mod_int(6)  = 0 ;

       intint_hot = fun_Gaussian_2D_double_integral_v1(x_int,y_int,par_mod_int);

       %- Summarize
       summary_quant_restrict(iD,1:4)  = data_det(iD,:);
       summary_quant_restrict(iD,5)    = intint_hot;
       summary_quant_restrict(iD,6)    = fit_spot.sigmaX;
       summary_quant_restrict(iD,7)    = fit_spot.amp;
       summary_quant_restrict(iD,8)    = fit_spot.bgd;
       summary_quant_restrict(iD,9)    = max(img_crop(:));
     end
end
disp(' .... DONE!')

%% Save results - free

%- File-name
[dum base] = fileparts(img_prot.file_names.raw);
file_name_default = fullfile(img_prot.path_names.img,[base,'__2D_QUANT_FREE.csv']);

% [file_save,path_save] = uiputfile(file_name_default,'Save results of HOTspots analysis');
% if file_save == 0; return; end  
% file_name_full = fullfile(path_save,file_save);
file_name_full = file_name_default;

%- Prepare data
N_par        = size(summary_quant,2);
cell_data    = num2cell(summary_quant);  
cell_write   = [cell_data];      
cell_write   = cell_write';   %- fprintf works on colums - data has to therefore be transformed 
string_write = ['%g',repmat(';%g',1,N_par-1),'\n'];

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'x;y;z;t;IntInt;Sigma;amp;bgd;IntMax\n');
fprintf(fid,string_write,cell_write{:});             
fclose(fid); 

%% Save results - restricted

%- File-name
[dum base] = fileparts(img_prot.file_names.raw);
file_name_default = fullfile(img_prot.path_names.img,[base,'__2D_QUANT_RESTRICT.csv']);

%[file_save,path_save] = uiputfile(file_name_default,'Save results of HOTspots analysis');
%if file_save == 0; return; end  
%file_name_full = fullfile(path_save,file_save);
file_name_full = file_name_default;

%- Prepare data
N_par        = size(summary_quant_restrict,2);
cell_data    = num2cell(summary_quant_restrict);  
cell_write   = [cell_data];      
cell_write   = cell_write';   %- fprintf works on colums - data has to therefore be transformed 
string_write = ['%g',repmat(';%g',1,N_par-1),'\n'];

%- write to file
fid = fopen(file_name_full,'w');
fprintf(fid,'x;y;z;t;IntInt;Sigma;amp;bgd;IntMax\n');
fprintf(fid,string_write,cell_write{:});             
fclose(fid); 