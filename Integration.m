% The model has three module's  - Sensitivity Module, Sky Module & Integration Module
% This is the third part of the model.
%%
% This is the called  Integration Module
% The input's to the module is Sensitivity Data and Sky Data obtained from simualtion of the previous two modules
% The output from the module is Irradiace on each cell of the PV module for every hour/one year and its corresponding PV module map
%% Module Description
% Sensitivity (LogFileRead) data and Sky data results are integrated here. 
% So run/load the Sensitivty data and Sky data before running this module.
% This module does not need any user input but only loading the sensitivity
% data and sky data.

%%
%INPUT Variables
%sky        sky (hemisphere)-radiance distribution form Sky model(normalized to horizontal diffuse irradiance & Direct Noraml Irradiance) [W/m2 SR]
%sens       surface (e.g. of PV system) sensitivity map showing power (or fraction of incident) absorbed for every incident direction divided by the irradiance of the emitting light source to give 'projected area' [no unit]

% OUTPUT Variables
% Cell_Irr   it is the information on the irradiance of the cells in the PV module for one day/year.  

%%
%Important:
    %Again "run the Sensitivity Module & Sky Module" or "load the Sensitivity Data & Sky Data" before running this module
%%
clearvars -except sen_interpolate str_sen Lv Time str1 str2 P


%for more than one cell
for a= 1:1:size(Lv,3)  % looping the no of hours in a day. If it is for one year then j=1, since size (Lv,3)=1
    Lv1=Lv(:,:,a);
for b = 1:1:size(sen_interpolate,3);          % no of cells in the module   - in this case it is 22 cells
       
sen1=sen_interpolate(:,:,b);


rp = size(Lv1,1) - 1;        %angular resolution (nr. intervals) in phi direction [-]
rz = size(Lv1,2) - 1;        %angular resolution (nr. intervals) in zeta direction [-]

zeta = linspace(0,pi/2,rz+1)'*ones(1,rp+1);

%area of sky element is circumference * width / rt

zeta_c = (zeta(1:end-1,1) + zeta(2:end,1)) / 2;    %center angle of zeta band
circumf = 2 * pi * sin(zeta_c) * ones(1,rp);        %circumference of zeta band (ones to make matrix)
width = pi / (2 * rz);                              %width of zeta band
dA = circumf * width / rp;                          %area of sky element

A = sum(sum(dA));                                   %total hemisphere area (should be 2*pi)
dA = dA * 2 * pi / A;
dA = dA';

sky_c = (Lv1(1:end-1,1:end-1) + Lv1(2:end,1:end-1) + Lv1(1:end-1,2:end) + Lv1(2:end,2:end)) / 4;
sens_c = (sen1(1:end-1,1:end-1) + sen1(2:end,1:end-1) + sen1(1:end-1,2:end) + sen1(2:end,2:end)) / 4;

%P11 = sky_c .* dA .* sens_c; 
P1 = sum(sum(sky_c .* dA .* sens_c));

P(b,a) = P1;    % the P matrix size here is 'b x a'
% Important Read:
% The columns of 'P matrix' correspond to the number of hours in a day. But for the year there will be only one column for P.
%The rows of 'P matrix' correspond to the number of cells. Here its is 22.
end
end

%P=round(P);   % roundiing of the matrix P
%climit = (max(max(P))-rem(max(max(P)),100))+100;   % to find the maximum irradiance of P for that particular day/year. This wiil be used in plotting.
%% Making the text appear in the plots
%{
% breaking str_sen to extract the strings to place them in the plots
str_break = strsplit(str_sen,'_');  

% Seperaion of the strings
str_break_mat_char = char(str_break(1));   % Material
str_break_mat_char_uppercase = upper(str_break_mat_char(1));
Material = [str_break_mat_char_uppercase,str_break_mat_char(2:end)];
%Material = strcat(cellstr(str_break_char_uppercase),cellstr(str_break_material_char(2:end)))
str_concatenated_material = (strcat('Reflector Materail = ',Material));
lth_concatenated_material = length(str_concatenated_material);

str_break_rang_char = (char(str_break(2)));  % Reflector Angle
str_rang = (str_break_rang_char(2:end));
str_concatenated_rang = (strcat('Reflector Angle = ',str_rang)); 
lth_concatenated_rang = length(str_concatenated_rang);

str_break_pvang_char = (char(str_break(3)));  % PV Angle
str_pvang = (str_break_pvang_char(3:end));
str_concatenated_pvang = (strcat('PV Angle = ',str_pvang)); 
lth_concatenated_pvang = length(str_concatenated_pvang);

str_break_ore_char = (char(str_break(4))); % Orientation
str_ore = (str_break_ore_char);
str_concatenated_ore = strcat('Orientation = ',str_break_ore_char);
lth_concatenated_ore = length(str_concatenated_ore);

max_lth_concatenated = max([lth_concatenated_material,lth_concatenated_rang,lth_concatenated_pvang,lth_concatenated_ore]);
str_concatenated_material = [str_concatenated_material,blanks(max_lth_concatenated - lth_concatenated_material)];
str_concatenated_rang = [str_concatenated_rang,blanks(max_lth_concatenated - lth_concatenated_rang)];
str_concatenated_pvang = [str_concatenated_pvang,blanks(max_lth_concatenated - lth_concatenated_pvang)];
str_concatenated_ore = [str_concatenated_ore,blanks(max_lth_concatenated - lth_concatenated_ore)];

%creating a chatacter cell array to be used in the placing of the text in the plot.
Plot_Text = [str_concatenated_material ; str_concatenated_rang ; str_concatenated_pvang ; str_concatenated_ore]
%}

%% 
curves = 1; % no of sine curves along the x-axis; one sine curve is 0 to 360 degrees - [no unit]
int = 5;
x  = ((0:int:360*curves) *pi/180)';   % radian points for creating a sine curve
z  = 0*sin(x);  
x1 = (x(1:end-1)+x(2:end))/2;
z1 = (z(1:end-1)+z(2:end))/2;
plot(x1,z1);
imagesc(x1,z1,P(:,2));

%%
%{
for a=2%1:1:size(Lv,3)
    l=1;
    m=1;
    
    PM_Time=P(:,a);
          
     for k=1:1:11
        for j=1:1:2
            Cell_Irr(j,k,a) = PM_Time(l,m,:);   % the data is in W/m2 .  The size of the 'Cell_Irr' data is '2x 11 x a'
            l=l+1;
        end
     end
     
    %Plotting
    figure(a+20)
    h(a)=imagesc(Cell_Irr(:,:,a));
    axis image     
        %Data inside the PV cell in the plot
    textStrings = num2str(Cell_Irr(:));  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    textStrings = textStrings(((22*a)-21):(22*a));
    [x,y] = meshgrid(1:(11),1:2);   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center','FontSize',14);   %# Plot the strings
    
    ax = gca; % current axes
    c = colorbar; % colorbar
    
    colormap hot
    caxis([0 1000])
    
    if strcmp(str2,'Full Year') == 1;     % for one year
        Time = [];
        ax.Title.String = strcat({'Irradiance on the PV Module for a Full Year in '},str1);
        %title(str,'FontSize',16)
        c.Label.String = 'kWh/m^2';
        else                                  % for one day
        ax.Title.String = strcat({'Irradiance on the PV Module on a '},str2,{' '},'in',{' '},str1,' at',{' '},num2str(Time(a)),' hrs');
        %title(str,'FontSize',16)
        c.Label.String = 'W/m^2';
    end
        % Making the text appear in the plots
        ylim=get(gca,'ylim');
        xlim=get(gca,'xlim');
        text(xlim(2)-11.10,ylim(1)-.5,Plot_Text,'color','r')   % here Plot_Text is a string matrix
   
    % assigining Axis lables for the plot
      ax.XTick = [1:11];
      ax.XTickLabel = {'1','2','3','4','5','6','7','8','9','10','11'};
      ax.YTick = [1:2];
      ax.YTickLabel = {'1','2'};
      ax.TickLength = [0 0];
      
    %% 
    %Maximizing plot
    %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
     
end

% saving the Data in a variable for further processing of that data by plotting various graphs.
%save(strcat('C:\Users\muthukumarva\Dropbox\3)SET\SEAC\Thesis\Lighttools-MATLAB\7) Final Model\Different Scenarios\4)Intergration Results\',str2,'\Int_',str_sen,'_',str1),'Cell_Irr','Time','str_sen','str1','str2');

%close all
%}
