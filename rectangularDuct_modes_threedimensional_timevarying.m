% This script produces an animation of the propagation of higher-order
% modes (as well as the plane wave) in a rectangular duct. It stores the
% result in a .avi file. 

clear variables

% length in three dimensions
Lx = 1;
Ly = 0.2;
Lz = 0.1;

% number of subdivisions along each dimension
Nx = 200;
Ny = 30;
Nz = 20;

% higher-order mode identification
m = 0;
n = 0;

% set alpha = 0 for no attenuation, vary it for increasing attenuation rate
% as in exp(-alpha*x) to simulate cut-off modes.

alpha = 0;

% normalizing factor for shape functions
    if (m == 0) && (n > 0)
        Nmn = sqrt(2);
    elseif (m > 0) && (n == 0)
        Nmn = sqrt(2);
    elseif (m > 0) && (n > 0)
        Nmn = 2;
    else
        Nmn = 1;
    end

kx = 10;
kz = n*pi/Lz;
ky = m*pi/Ly;

% common to all faces
x = linspace(0,Lx,Nx);

% Specify initial time (t_0), final time (t_end) and number of subdivisions in time (N_time) 
% for the animation

t_0 = 0;
t_end = 1;
N_time = 100;
time = linspace(t_0,t_end,N_time);
loops = numel(time);

% specify a frequency in [Hz]
f = 2;
omega = 2*pi*f;

% name of the movie where the animation is stored
v = VideoWriter('00_rectangle_duct.avi');
v.FrameRate = 30;

open(v)

for indx = 1:loops
    
    % face #1 (y = 0)

    z = linspace(0,Lz,Nz);

    [X,Z] = meshgrid(x,z);
    
    % create an interconnection matrix between the points in the grid
    T = delaunay(X,Z);

    Y = zeros(size(X));
    
    % create a matrix for the coordinates of the vertices

    V = [X Y Z];

    % Rewrite it so that each row corresponds to the coordinates (x,y,z) of
    % the i-th vertex
    V = reshape(V,numel(X),3);

    % Compute the shape function in each vertex
    C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1)).*exp(-alpha*V(:,1));


    figure(1)

    % create a patch specifying the faces and the vertices, as well as the
    % value of the computed function at the vertices
    h = patch('faces',T,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp');
    set(h,'FaceAlpha',1,'EdgeColor','none')
    
    hold on
    
    % repeat for all the faces

    % face #2 (z = 1)
    y = linspace(0,Ly,Ny);

    [X,Y] = meshgrid(x,y);
    T = delaunay(X,Y);

    Z = Lz*ones(size(Y));

    V = [X Y Z];

    V = reshape(V,numel(Y),3);

    C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1)).*exp(-alpha*V(:,1));

    h = patch('faces',T,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp');
    set(h,'FaceAlpha',1,'EdgeColor','none')

    % face #3 (y = Ly)

    z = linspace(0,Lz,Nz);

    [X,Z] = meshgrid(x,z);
    T = delaunay(X,Z);

    Y = Ly*ones(size(X));

    V = [X Y Z];

    V = reshape(V,numel(X),3);

    C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1)).*exp(-alpha*V(:,1));

    h = patch('faces',T,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp');
    set(h,'FaceAlpha',1,'EdgeColor','none')

    % face #4 (z = 0)
    y = linspace(0,Ly,Ny);

    [X,Y] = meshgrid(x,y);
    T = delaunay(X,Y);

    Z = zeros(size(Y));

    V = [X Y Z];

    V = reshape(V,numel(Y),3);

    C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1)).*exp(-alpha*V(:,1));

    h = patch('faces',T,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp');
    set(h,'FaceAlpha',1,'EdgeColor','none')
   
    % face #5 (x = Lx)
    y = linspace(0,Ly,Ny);
    z = linspace(0,Lz,Nz);

    [Y,Z] = meshgrid(y,z);
    T = delaunay(Y,Z);

    X = Lx*ones(size(Y));

    V = [X Y Z];

    V = reshape(V,numel(X),3);

%     C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1));
    C = Nmn*cos(kz*V(:,3)).*cos(ky*V(:,2)).*cos(omega*time(indx)-kx*V(:,1)).*exp(-alpha*V(:,1));

    h = patch('faces',T,'Vertices',V,'FaceVertexCData',C,'FaceColor','interp');
    set(h,'FaceAlpha',1,'EdgeColor','none')
    
    % Set dimension of the image
    
    set(gcf,'Units','Inches');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','Position',[1,1,15,15]);
    
    % custom colorbar

    J = customcolormap_preset('red-white-blue');
    colormap(J);

    axis off

    box("off")
    
    % set aspect ratio

    pbaspect([Lx/Lx Ly/Lx Lz/Lx]);

    % set viewing position and angles
    caxis([-1.5 1.5]);
    camtarget([Lx/2,Ly/2,Lz/2])
    campos([4,-3,2.5])
    camup([0,0,1]);
    camva(10);
    view([45, 25]);
    
    % get the frame and store it in a video
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    close all
    
end

close(v)