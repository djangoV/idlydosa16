% The model has three block's  - Sensitivity Block, Sky Block & Integration Block
% This is the second block of the model.

%%
% This is called the Sky Block
% The input's to this module is the location's Azimuth & Zenith angles of the sun and its corresponding DHI & DNI
% The output's from this module will be Sky radiance distribution or Sky Data(variable is "Sky_data" and Sky_Time) and the corresponding Sky maps. 

%% Block Description
% South Azimuth = 0 degrees
% The range of the angles --> Azimuth = -180 to +180 & Zenith = 90 to 0
% The angle step used here is 0.5 degrees i.e) Full Sky Scan
% For the Direct component the sun is placed in the sky based on its position
% For the Diffuse component the Perez Luminance Model is used
% Direct + Diffuse will give the total radiance at that specific point in the sky.

%% Important Read - Description for Inputing data in a specified format
% The Data is obtained from Meteonorm - year 2005
% It should only be specified in the below mentioned format
% The data format type is in excel and it contains the following columns "year,month,day,hour,Azimuth,Altitude,DNI,DHI,GHI,Temperature".
% Also the perez coefficients should be loaded in to the same file as the simulation file. - Perez_table.mat

%% Note 1:
% 1.1)
        % "int"  = Interval step for the sky distribution. 
        % 0.5 degrees is the only choice because of the accurate prediction of the sun in the sky for an hour.
        % This is because the sun's angular diameter is 0.5 degrees and hence the sky discretization interval is also 0.5 degrees.

% 1.2) Abbrevations used in this module:
    % "Az"  = Azimuth
    % "Z"   = Zenith
    % "Ees" = Direct Normal Irradiance
    % "Eed" = Diffuse Horizontal Irradiance

% 1.3) Important - This module is used for simulation of the Sky for every hour in one day of the year.  

%% Note 2:
% 2.1)You can either run the Sensitivity Module before or after this module

%% The block starting:
clearvars 

% Note 3:
% 3.1)In some cases the Lv_diffuse produces complex data, hence only real part of the data is used. 

%% Fixed Inputs
months = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'June' 'July' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};  % this will be used in naming variables

%% Inputting excel File
str1 = num2str('Phoenix');  % User Input - City
ExcelFile = strcat('C:\Users\Adi\Dropbox\3)SET\SIP2\Lux\LUX_2016_04_14\Sky Data\Meteonorm-',str1,'-year.xlsx');
ExcelData=xlsread(ExcelFile);

%% Filtering

for p = 3% 1:1:12;        % no of months in a year
    
    ExcelData_filter_1month = ExcelData(ExcelData(:,2)== p,:);  % filering data for one month
    days = length(unique(ExcelData_filter_1month(:,3)));       % finding out how many days are there in the month.
        
    for q = 28%1:1:days       % no of days in a month
        
        field = char(strcat(months(p),'_',num2str(q)));            % used in structural array
        
        ExcelData_filter_1day = ExcelData_filter_1month(ExcelData_filter_1month(:,3)== q,:);  % filtering data for one day
        
        %Filtering the data for sun aboive horizon
        ExcelData_filter = ExcelData_filter_1day(ExcelData_filter_1day(:,6)>0,:);   % filering data for sun above the horizon i.e) sun altitude>0
        Z= 90-ExcelData_filter(:,6);                                                % corresponding sun zenith angle.   --> here the altitude is changed to zenith angle
        Az= ExcelData_filter(:,5);                                                  % sun azimuth --> input is South = 0 degrees

        % Sun Rise hour and Sun Set hour and the times in-between them
        time=min(ExcelData_filter(:,4)):1:max(ExcelData_filter(:,4));               % Time during which the sun is above horizon in the day (hrs).
        

        % Irradiance Input ; Both Diffuse and Direct        
        Ees0 = 1361;                                                                % normal incident extraterrstrial irradiance [W/m2]

        % DNI & DHI Input
        Ees= ExcelData_filter(:,7);                                                 % Direct Normal Irradiance.
        Eed= ExcelData_filter(:,8);                                                 % Diffuse Horizontal Irradianace

        % Dividing the points in the Sky Hemisphere
        int = 0.5;                                                                  % angle step - this is an input.

        %{
        % this is for reference -units in degrees and will not be used in the module
        t0 = (-180:int:180)';                                                       % sky grid azimuth angle vector [deg ]
        z0= (0:int:90)';
        [theta0,zeta0] = meshgrid(t0,z0); 
        %}

        %this unit is in radians - will be used in the block
        t = (-180:int:180)'*pi/180;                                                 % sky grid azimuth angle vector [rad]
        z = (90:-int:0)'*pi/180;                                                    % sky grid zenith angle vector [rad]
        [theta,zeta] = meshgrid(t,z);                                               % sky grid zenith and azimuth angle matrix [rad]

        % defining the grid in the hemisphere
        rt = size(zeta,2) - 1;        %angular resolution (nr. intervals) in theta direction [-]
        rz = size(zeta,1) - 1;        %angular resolution (nr. intervals) in zeta direction [-]

        % Length of the "Excel Data_filter" 
        lth = size(ExcelData_filter,1);                                               % calculate the length of the rows in the Excel Data_filter variable.
                                                                      

        %Creating the Empty Luminance Array
        Lv=zeros(rt+1,rz+1,lth);

        %Creating Empty Arrays for adding data from every hour
        Lv_Direct_Final =zeros(rz+1,rt+1); 
        Lv_Diffuse_Final =zeros(rz+1,rt+1);



        for i = 13%:1:lth              % no of sun-hours in a day
    
            Az1 = Az(i)*pi/180;             % Azimuth loaded from meteonorm data
            Z1 = Z(i)*pi/180;               % Zenith loaded from meteonorm data

            Ees1 = Ees(i); %DNI = W/m2    % Direct For One Hour from meteonorm data
            Eed1 = Eed(i); % DHI = W/m2;  % Diffuse For One Hour from meteonorm data


            %--------------------------------------------------------------------------%--------------------------------------------------------------------------   
            gamma = acos(cos(Z1) * cos(zeta) + sin(Z1) * sin(zeta) .* cos(Az1-theta));   % [rad]
            %--------------------------------------------------------------------------%-------------------------------------------------------------------------- 

            %---sky clearness---
            epsilon = ((Eed1 + Ees1)/Eed1 + 1.041*Z1^3) / (1+1.041*Z1^3);  %
            %eq. (4) of ref.1

            %---sky brightness---
            m = 1/cos(Z1);                  %optical air mass [-] 
            delta = m * Eed1 / Ees0;  

            %--------------------------------------------------------------------------%--------------------------------------------------------------------------

            load('Perez_Table.mat','A','B','C','D','E'); %load table from file

            % check for the irradainces  - value to be positive, if negative then it should be zero.

            %---select appropriate table row based on sky clearness---
            if epsilon <1.000
                error('sky clearness < 1, outside range of Perez table 1.');
                %should never happen
            elseif epsilon < 1.065
                row = 1;
            elseif epsilon < 1.230
                row = 2;
            elseif epsilon < 1.500
                row = 3;
            elseif epsilon < 1.950
                row = 4;
            elseif epsilon < 2.800
                row = 5;
            elseif epsilon < 4.500
                row = 6;
            elseif epsilon < 6.200
                row = 7;
            else
                row = 8;
            end

            %--- empirical relations a = f(a1,a2,a3,a4), same for b,c,d,e---
            %subcoefficients are taken from loaded data
            a = A(row,1) + A(row,2)*Z1 + delta * (A(row,3) + A(row,4)*Z1);         %eq(6)
            b = B(row,1) + B(row,2)*Z1 + delta * (B(row,3) + B(row,4)*Z1);         %eq(6)
            e = E(row,1) + E(row,2)*Z1 + delta * (E(row,3) + E(row,4)*Z1);         %eq(6)

            if row == 1
                c = exp((delta * (C(row,1) + C(row,2) * Z1))^C(row,3)) - C(row,4);    %eq(7)
                d = -exp(delta*(D(row,1)+D(row,2)*Z1)) + D(row,3) + delta * D(row,4); %eq(8)
            else
                c = C(row,1) + C(row,2)*Z1 + delta * (C(row,3) + C(row,4)*Z1);         %eq(6)
                d = D(row,1) + D(row,2)*Z1 + delta * (D(row,3) + D(row,4)*Z1);         %eq(6)
            end
            %--------------------------------------------------------------------------%--------------------------------------------------------------------------       
            exp_cos_zeta = exp(b ./ cos(zeta));
            lv_diffuse = (1 + a * exp(b ./ cos(zeta))) .* (1 + c * exp(d * gamma) + e * cos(gamma).^2); %eq(1)
            % lv_diffuse(isinf(lv_diffuse)) = 0;
            %--------------------------------------------------------------------------%--------------------------------------------------------------------------

            % Diffuse Irradiance Normalization

            %sky element central value is average of four corner grid points
            lv_c_diffuse = (lv_diffuse(1:end-1,1:end-1) + lv_diffuse(2:end,1:end-1) + lv_diffuse(1:end-1,2:end) + lv_diffuse(2:end,2:end)) / 4;

            %area of sky element is circumference * width / rt


            zeta_c = (zeta(1:end-1,1) + zeta(2:end,1)) / 2;     %center angle of zeta band
            circumf = 2 * pi * sin(zeta_c) * ones(1,rt);        %circumference of zeta band (ones to make matrix) when u visualize in a 2D polar chart.
            width = pi / (2 * rz);                              %width of zeta band
            dA = circumf * width / rt;                          %area of sky element

            A = sum(sum(dA));                                   %total hemisphere area (should be 2*pi)
            dA = dA * 2 * pi / A;                               %minor correction (due to discretization) % nomalization

            I_diffuse = sum(sum(lv_c_diffuse .* dA .* (cos(zeta_c) * ones(1,rt))));  %irradiance when not normalized
            %note: cos(zeta) is included because *horizontal* diffuse irradance 
            %is proportional to lv * cos(zeta)! (this explains eq(3) in ref.1)

            %normalize such that irradiance becomes equal to Eed
            Lv_Diffuse = lv_diffuse * Eed1 / I_diffuse;

            Lv_Diffuse = flipud(Lv_Diffuse);
            Lv_Diffuse = real(Lv_Diffuse);     % to remove any complex part if arises. Yes they do aries for monthly simulation
            Lv_Diffuse_Final = Lv_Diffuse_Final+Lv_Diffuse;

            %--------------------------------------------------------------------------%--------------------------------------------------------------------------

            % Direct Irradiacne and its Normalization
            lv_direct = zeros(length(z),length(t));

            Az_1 = Az(i) + 180; 

            Var1=(round(Az_1/int))+1;     % Azimuth Angle
            Var2=(round(Z(i)/int))+1;     % Zenith Angle

            lv_direct(Var2,Var1) = 1;

            lv_c_direct= (lv_direct(1:end-1,1:end-1) + lv_direct(2:end,1:end-1) + lv_direct(1:end-1,2:end) + lv_direct(2:end,2:end)) / 4;

            I_direct = sum(sum(lv_c_direct.* dA));  %irradiance when not normalized

            Lv_Direct = lv_direct * Ees1 / I_direct;   % Normalizing the direct Irradiance

            Lv_Direct_Final = Lv_Direct_Final + Lv_Direct;  % 90 to 0 zenith

            Lv1 = Lv_Direct + Lv_Diffuse; 
            Lv1=flipud(Lv1);   % flipping from 90 to 0 zenith, to ,0 to 90 altitude. 
            Lv1=Lv1';   % matix transpose from 181x721 to 721x181
            Lv(:,:,i) = Lv1(:,:);    % Final Variable that will be used in the Integration Module

            Sky_data.(field) = Lv;    % The Sky Data output variable stored as a structure array.
            Time_Sky.(field) = time;    % The Time of the day when the sun is above horizon. It is also a structure array.
        end
    end
end

% Suming up the radiances for a certain period (say day, month,year)

names = fieldnames(Sky_data);
structSize = length(names);

%Creating Empty Arrays for adding data from every hour
Lv_Final =zeros(rt+1,rz+1); 

for p = 1:1:structSize   
Sky_data1 = Sky_data.(char(names(p)));
Lv_Final = Lv_Final+sum(Sky_data1,3);
end
Lv_Final = Lv_Final./1000;   % to divide it by 1000 to covert Wh/m^2/sr to kWh/m^2/sr.

%%
%{
%Plotting
% conversion into spherical coordinates to be used in the Dome plots.
XX = cos(zeta) .* cos(theta);
YY = cos(zeta) .* sin(theta);
ZZ = sin (zeta) ;

% conversion into polar coordinates to be used in the polar plots.
X = zeta .* sin(theta);
Y = zeta .* cos(theta);

%plots

%polar plot - for 2D View -Direct + Diffuse 
figure(98)
pcolor(X,Y,Lv_Final')
set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
axis equal tight off
shading flat
text (0,1.68,0,'S','FontSize', 30)
text (0,-1.68,0,'N','FontSize', 30)
text (1.8,0,0,'W','FontSize', 30)
text (-1.555,0,0,'E','FontSize', 30)
cc = colorbar;
cc.Label.String = 'kWh/m^2/sr';
cc.Label.FontSize = 45;
cc.FontSize = 30;
caxis([0 400])


%}
