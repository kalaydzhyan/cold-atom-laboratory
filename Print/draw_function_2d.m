%% Draw wave function
%% INPUTS:
%%          phi: Wave function (matrix)
%%          Method: Structure containing variables concerning the method (structure) (see Method_Var2d.m)
%%          Geometry2D: Structure containing variables concerning the geometry of the problem in 2D (structure) (see Geometry2D_Var2d.m)
%%          Figure: Structure containing variables concerning the figures (structure) (see Figure_Var2d.m)

function draw_function_2d(phi,Geometry2D,Figure)

global glTime;
global gldeltaT;
global Timeline;
global now_stamp;
global color_limit;

figure(Figure.label); % Setting the number of the figure
clf(Figure.label); % Clear figure
pcolor(Geometry2D.X,Geometry2D.Y,phi); % Drawing function

if length(color_limit) == 0
    caxis([0 1.3E-3]);
else
    caxis([0 color_limit]);
end

%contourf('v6',Geometry2D.X,Geometry2D.Y,phi); % Drawing function
axis equal
axis tight
set(gcf,'color','w');
shading interp; % Setting shading
colormap(Figure.map); % Setting colormap
colorbar; % Setting colorbar
view(2); % Setting view
%xlabel(text, 'FontSize', 30) 
xlabel(Figure.x, 'FontSize', 17); % Setting x-axis label
ylabel(Figure.y, 'FontSize', 17); % Setting y-axis label
title(Figure.title, 'FontSize', 17); % Setting title of the figure
set(gca,'fontsize',17)
if (isvector(Figure.axis)==1)
    caxis(Figure.axis)
end

print(strcat(now_stamp,'/', num2str(glTime, 5),'.png'), '-dpng');
Timeline = [Timeline glTime];
glTime = glTime + gldeltaT;
drawnow; % Drawing
