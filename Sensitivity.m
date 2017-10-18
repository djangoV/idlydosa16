
% clearvars
 %% Hyat Module Dimensions: Not used in the model.
%{
% It is given here only for Contex purpose
%Length of undulated module is in m;
length = 5.91;
%Width of the undulated module is in m;
width = 0.3;
%period of the undulated module is in m;
period = 0.19;
%amplitude of the undulated module is in m;
amplitude = 0.02
%}
%% User Input Variables 
% PV/ BIPV Undulation Design:
lth = 0.19;                             % length of the curve - [meters]
amp = 0.02;                                % amplitude of the sine wave - [meters]
prd = 0.19;                             % length of period of the module [meters]
x_pts = 51;                             % the number of points to discretize the sine curve along the x-axis - [no unit]        
wth = 0.1;                              % width of the module - [meters]
y_boxes = 1;                            % the number of boxes along the y-axis - [no unit]

%Orientation
Direction = 0;  % this is the input by the user.

% Hemisphere Scan
az_int = 10;                            % azimuth interval to perform the dome scan - [degrees]

%% Fixed Input Dependant Variables
% PV/ BIPV Module Design: - For Sine curve that will used in the script

curves = lth/prd;                       % no of periods/undulations in the sine curve - [no unit]           
wth_boxes  = wth/y_boxes;               % width of the boxes [meter]
lambda = (prd/pi)*pi;                   % period of the sine curve  - [radians]
omega = (2*pi)/(lambda);                % [no unit]
X_limit = (360*curves)/omega;           % [degrees]
int = X_limit/(x_pts- 1);               % [degrees]
Y_wth = wth_boxes* y_boxes;             % [meter]

%% Creation of the radian points (x and z coordinates) of the sine curve

x  = (0:int:X_limit) *pi/180;                   % x coordiantes for creating the sine curve - [radians]
z  =  amp*sin(omega*x);                         % z coordinates for creating the sine curve - [radians]
L_curve = sum(sqrt( diff(z).^2 + diff(x).^2));  % Calculation of length of the sine curve

L_seg = L_curve/(size(x,2)-1);                  % Length of each arc in the curve

% L_arc = prd/(size(x,2)-1);

% % Testing 1
%  figure(100)
%  plot(x,z)     % this plot is plotted to view the look of the since curve 

%% Equal Arc length
% This section is about calculating the new x coordinates that makes the arc lengths of the sine curve equal in length. 


x_new = zeros(1,size(x,2));                                                 % pre allocation of the new x coordinates. It is the same size as the variable 'x'                 


for a = 1:1:(size(x,2)-1)                                          

    l_ite=0;                                                                % pre allocation of the 'l_ite' variable.                                                
    x1 = x_new(a);
    
    while l_ite < L_seg                                                     % iteration is used to generate a new x coordiante.

        x1 = x1 + ((X_limit*pi)/180)/((size(x,2)-1)*10000);                 % calculation of the new x coordinate by iteration. 
                                                                            % 10000 in the denominatior is used to increase the calculation accuracy.
                                                                            
        l_ite  = sqrt((amp*sin(omega*x1)-amp*sin(omega*x_new(a)))^2+(x1-x_new(a))^2);   % length calculation. 
    
    end

         x_new(a+1) = x1;           % the new x coordinates are assigned to the variable 'x_new'
    
        % % Testing 2
        % L_arc_new(a) = l_ite;     %    the new length (equal) of each arc in the curve is stored in this variabe.
   
end

z_new  =  amp*sin(omega*x_new);     % the new z coordinates for creating the sine curve, after making all the arc lengths equla in size.

% % Testing 3
% L_curve_new = sum(sqrt( diff(z_new).^2 + diff(x_new).^2));  % Calculation of length of the new sine curve
% % Result  - Same as the original sine curve.

% % Testing 4
% figure(99)
% plot(x_new,z_new);  % This plot is created to view the new sine cuve that has equal arc lengths 


%% % Coordinates For the Lux Software. 
%V = cat(1, zeros(8,3), zeros((y_boxes+1)*size(x,2),3));   % Pre allocating the size of V array, ie the number of points to be plotted.
%F = cat(1, zeros(6,4),zeros (y_boxes*(size(x,2)-1),4));   % Pre allocating the size of F array, i.e the number of facets to be created.

%%
y1 = zeros(1,size(x,2));            % the y coordinates (width) of the sine curve. The width of the sine curve starts from y axis = 0;

V1 = cat(1,x_new,y1,z_new)';        % the xyz coordiantes of the vertices for the first since curve seen from the fromt view (xz plane).

V_2(1:size(x,2), 1:3) = V1; 

 
for b = 1:1:y_boxes;
    y2 = b * wth_boxes*ones(1,size(x,2));      % the width of the boxes in y axis is for every wth [no unit] interval.
    V2 = flipud(V1);                   
    V2(:,2) = y2(:);
    V1 = V2;                        % data points created for the second, third(looped) ...  sine curves seen from the front view.
    V_2(((size(x,2)*b)+1):(b+1)*size(x,2),1:3) = V2;   % inserting the second, third,.. since surve data points to the final array V
    F_new = (cat(1,((b-1)*size(x,2)+1):((size(x,2)*b)-1), ((b-1)*size(x,2)+2):(size(x,2)*b) , ((size(x,2)*(b+1))-1):-1:((size(x,2)*b)+1) , (size(x,2)*(b+1)):-1:((size(x,2)*b)+2)))';
    F_2(((b-1)*(size(x,2)-1)+1):((size(x,2)-1)*b),1:4) = F_new; % creating facet data points 
end


 %Box vertices
 V_1 =[
      (X_limit)*pi/180  0        -(amp+0.01) ;               %1 data point
      (X_limit)*pi/180  0        (amp+0.01)  ;               %2      "
      (X_limit)*pi/180  Y_wth    (amp+0.01)  ;               %3      "
      (X_limit)*pi/180  Y_wth    -(amp+0.01) ;               %4      "
      0                 Y_wth    -(amp+0.01) ;               %5      "
      0                 Y_wth    (amp+0.01)  ;               %6      "
      0                 0        (amp+0.01)  ;               %7      "
      0                 0        -(amp+0.01)];               %8      "



% Box facets
F_1 = [ 6  7  2  3;   %1 box - ceiling
        5  8  1  4;   %2 box - floor 
        1  2  3  4;   %3 box - side walls
        3  4  5  6;   %4 box - side walls
        6  7  8  5;   %5 box - side walls
        7  2  1  8 ]; %6 box - side walls
  
 %% simuluation run i.e) creation of the points and facets
 V = cat(1, V_1, V_2);
 F = cat(1, F_1, (F_2+8));  % addition of +8 because of the eight data points created fot the box around the sine curves.
%%
Lux(V);
Lux(V,F);
%% Assigning surface optical properties and teleportaion boundaries

Type(1).RT    = {0 0 'floor'};         %absorbing floor
Type(1).Plot  = {[0.3 0.3 0.3] 0.3};   %gray
Type(1).Facet = 2;

Type(2).RT    = {0 1 'wall'};          %transmitting side walls
Type(2).Plot  = {[0.0 0.0 1.0] 0.1};   %blue
Type(2).Facet = 3;
Type(2).Teleport = lambda;                      

Type(3).RT    = {0 1 'wall'};          %transmitting side walls
Type(3).Plot  = {[0.0 0.0 1.0] 0.1};   %blue
Type(3).Facet = 4;
Type(3).Teleport = -Y_wth;

Type(4).RT    = {0 1 'wall'};          %transmitting side walls
Type(4).Plot  = {[0.0 0.0 1.0] 0.1};   %blue
Type(4).Facet = 5;
Type(4).Teleport = lambda;

Type(5).RT    = {0 1 'wall'};          %transmitting side walls
Type(5).Plot  = {[0.0 0.0 1.0] 0.1};   %blue
Type(5).Facet = 6;
Type(5).Teleport = Y_wth;

Type(6).RT    = {0.12 0 'sinewave#'};     %absorbing surface
Type(6).Plot  = {[1 0 0] 0.7};         %red
Type(6).Facet = [7:size(F,1)];

%% Sensitivity Block 

% Hemisphere Scan 
azi = -180:az_int:180;                          % azimuth of the light source during the dome scan
ele = [89.5 90-[10:10:90]];                     % elevations/altitudes of the light source during the dome scan

sen = zeros(size(azi,2),size(ele,2),size(F_2,1)); % preallocating the size of sensitivity data array to reduce the speed

% Area of the Light Source and Cell:
A_LS = ((X_limit)*pi/180) * Y_wth;
A_seg = L_seg * wth_boxes;
area_ratio = A_LS/ A_seg;


boxes_even_no = floor((y_boxes/2))*2;            % (number of even-number boxes) * 2 


for i = 1:1: size(azi,2)                                % azimuth angle  South  = -ve x-axis.
    for j = 1:1:size(ele,2)                             % altitude angle 
        Type(7).RT    = {0 0 'ceiling'};                % absorbing ceiling
        Type(7).Plot  = {[0.7 0.7 0.7] 0.3};            % light gray
        Type(7).Emit  = {62500 [ele(j) ele(j) ; azi(i) azi(i)] 0};% {ceiling no of rays  [altitude altitude ; azimuth azimuth] random number } South  = -ve x-axis.
        Type(7).Facet = 1;

        %% Simulation Run for one source light position
        wl = 0.600;             %600 nm wavelength
        [OUT,T] = Lux(V,F,Type,wl);
        % pause(2)
        % view([90,0])

        lux_out = OUT(6).facc(:,2);

        % to flip the lux_out data of the even numbered cells in the module.
        if boxes_even_no ~= 0
            for k = 2:2:boxes_even_no
                lux_out(((size(lux_out,1)/y_boxes)*(k-1)+1):((size(lux_out,1)/y_boxes)*k)) = flipud(lux_out(((size(lux_out,1)/y_boxes)*(k-1)+1):((size(lux_out,1)/y_boxes)*k)));   % flipping the data of second, fourth, sixth,... cell strips in the undulated module.
            end
        else
            lux_out=lux_out;
        end
        sen_out = permute(lux_out, [3 2 1]).* area_ratio .* sind(90-ele(j));  % rearranging the LUX OUTPUT data as an 3D matrix 1x1xsize(facets) and correcting the Lux output by including the sine effect and area of the light source, cell area. 
        sen(i,j,:) = sen_out;                                                 % sensitivity data of the surfaces after the ray tracing simulation size(i)xsize(j)xsize(facets)

    end
end


%% interpolation of the sensitivity data
ze = (0:10:90);                          %sky grid altitude angle vector [rad]
t = (-180:az_int:180);                  %sky grid azimuth angle vector [rad]
[theta,zeta] = meshgrid(t,ze);
%%
zeze = (0:0.5:90);                        %sky grid altitude angle vector [rad]
tt = (-180:0.5:180);                    %sky grid azimuth angle vector [rad]
[theta_1,zeta_1] = meshgrid(tt,zeze);
%%
for s = 1:1:size(sen,3) 
sen_interpolate(:,:,s) = griddata(theta',zeta',sen(:,:,s),theta_1',zeta_1');
end
%%
%load('h:\Desktop\SIP-2 - Dropboxed 29-7-2016\Lux\LUX_2016_04_14\Simulation Run\6)Hyet_sinx_surface_30-7-2016\sen_interpolate_hyet__surface.mat');

%% Orientation sub-block
% Direction = 45;

if Direction == 0;
    Index = 0;
    Sensitivity = circshift(sen_interpolate,Index);
else
Index = 2*Direction; 
if Index > 0
    Index = Index +1;
else
    Index = Index;
end
Sensitivity = circshift(sen_interpolate,Index);
end
%}
%% Plotting

%----------------------------------------------------------------------------------------------------------------
% Before Interploation Interval - Polar Coordinates
aa = sind(t') * (-ze);
bb = cosd(t') * (-ze);
%%
% After Interpolation Interval - Polar Coordinates
aa_Interploate = sind(tt') * (-zeze);
bb_Interpolate = cosd(tt') * (-zeze);
%%
for s = 1%:1:size(sen,3)
% Polar Plot - before interpolation
%{
    figure(s) 
    pcolor(aa,bb,fliplr(Sensitivity(:,:,s)))
    axis equal tight off
    hold on
    text (0,94,'N')
    text (0,-94,'S')
    text (90.3,0,'E')
    text (-96.5,0,'W')
    shading interp
    colorbar
    colormap hot

% Polar Plot - after interpolation and before orientation
    figure(s+50)                         %new figure
    pcolor(aa_Interploate,bb_Interpolate,fliplr(Sensitivity(:,:,s)))
    axis equal tight off
    hold on
    text (0,94,'N')
    text (0,-94,'S')
    text (90.3,0,'E')
    text (-96.5,0,'W')
    shading interp
    colorbar
    colormap hot
%}
% Polar Plot - after interpolation and after orientation
    figure(s+100)                         %new figure
    pcolor(aa_Interploate,bb_Interpolate,fliplr(Sensitivity(:,:,s)))
    axis equal tight off
    hold on
    text (0,96,'N','FontSize',28)
    text (0,-96,'S','FontSize',28)
    text (91,0,'E','FontSize',28)
    text (-103,0,'W','FontSize',28)
    shading interp
    cc = colorbar;
    colormap hot
    cc.FontSize = 42;

end
%before_orientating_sen_interpolate = sen_interpolate;


%% Figure Background color
%{
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1])         % [0 0 0] - black ; [1 1 1] - white
%}
%}
%% Sky Block
% Sky_Block_Meteonorm;
