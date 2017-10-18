%% Sky Map plot
%% Input
% Load the "Lv_Final" variable obtained from the "Sky_Block_Meteonorm.m" MATLAB file.
% Here the "LV_Final" should be of array size 721x181, going from 0 to 90 altitude and -180 to +180 azimuth

%% Output
% It is the sky map.
% Note:
    % If the sky map is for an instant in time, the caxis unit is W/m^2/sr.
    % If the sky map is for an integrated period, the caxis unit is
    % Wh/m^2/sr or kWh/m^2/sr.

%% Script Plotting
    int = 0.5;   % interval of the sky dome.
%this unit is in radians - will be used in the block
    t = (-180:int:180)'*pi/180;                                                 % sky grid azimuth angle vector [rad]
    z = (90:-int:0)'*pi/180;                                                    % sky grid zenith angle vector [rad]
    [theta,zeta] = meshgrid(t,z);

% conversion into spherical coordinates to be used in the Dome plots.
    XX = cos(zeta) .* cos(theta);
    YY = cos(zeta) .* sin(theta);
    ZZ = sin (zeta) ;

% conversion into polar coordinates to be used in the polar plots.
    X = zeta .* sin(theta);
    Y = zeta .* cos(theta);

% plots
% polar plot - for 2D View -Direct + Diffuse 
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
    cc.Label.String = 'kWh/m^2/sr';   % always check the unit. see above in "%% output"
    cc.Label.FontSize = 45;
    cc.FontSize = 30;
    caxis([0 400])

% 3D Dome - Direct or Diffuse
    figure(99)
    surf (XX,YY,ZZ,flipud((Lv_Final')))
    axis equal tight
    shading interp
    colorbar
    colormap parula
    view(67,18)
    % view(90,90)
    grid off
    zlabel('Elevation','FontSize', 16)
    set(gca,'Ydir','reverse')
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    set(gca,'ZTickLabel','')
    text (0,1.25,0,'W','FontSize', 16)
    text (0,-1.25,0,'E','FontSize', 16)
    text (-1.25,0,0,'N','FontSize', 16)
    text (1.25,0,0,'S','FontSize', 16)
    caxis([0 400])          % always check the unit. see above in "%% output"
