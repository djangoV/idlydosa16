%Plotting
%%%%%%%%%%%%% SENSITIVITY MAP %%%%%%%%%%%
%Plotting
%----------------------------------------------------------------------------------------------------------------
zz = (0:0.5:90);          %sky grid zenith angle vector [rad]
tt = (-180:0.5:180);         %sky grid azimuth angle vector [rad]
[theta_1,zeta_1] = meshgrid(tt,zz);

% After Interpolation Interval - Polar Coordinates
aa_Interploate = sind(tt') * (-zz);
bb_Interpolate = cosd(tt') * (-zz);

% After Interpolation Interval - Spherical Coordinates Coordinates
[theta_1,zeta_1] = meshgrid(tt,zz);
xx_Interpolate = (cosd(zeta_1) .* cosd(theta_1))';
yy_Interpolate = (cosd(zeta_1) .* sind(theta_1))';
zz_Interpolate = (sind (zeta_1))';

%Plots
 for s = 11 % 1:1:22      %          s =11 denotes the center cell in the PV module.
     %for every simulation
      
 % After Interpolation- Polar Plot      
    figure(s+50)                         %new figure
    pcolor(aa_Interploate,bb_Interpolate,fliplr(sim_interpolate(:,:,s)))
%             xlabel('Altitude') % x-axis label
%             ylabel('Altitude')
%             title(str(s))
    axis equal tight off
    hold on
    text (0,95,'N','FontSize', 25)
    text (0,-95,'S','FontSize', 25)
    text (88.9,0,'E','FontSize', 25)
    text (-103,0,'W','FontSize', 25)
    shading interp
    %colorbar('FontSize',30);
    
    c = colorbar;
%     c.Label.String = '[-]/sr';
%     c.Label.FontSize = 18;
     c.FontSize = 14;
    colormap hot
  
 end

