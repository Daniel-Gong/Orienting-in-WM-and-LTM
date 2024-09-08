
function [colourCircle,colourBar] = create_colour_circle(fullcolormatrix)

format long g;
format compact;
fontSize = 24;
scaling = 0.7;
% width = 1440;
height = 900;
threesixty = 0; %switch to 0 for 180deg hsv
%==================================================================================================
%% Compute HSV image
rows = height-height*scaling;
columns = height-height*scaling;
midX = columns / 2;
midY = rows / 2;
% Construct v image as uniform.
v = 0.95 * ones(rows, columns);
s = zeros(size(v)); % Initialize.
h = zeros(size(v)); % Initialize.
% Construct the h image as going from 0 to 1 as the angle goes from 0 to 360.
% Construct the S image going from 0 at the center to 1 at the edge.
for c = 1 : columns
	for r = 1 : rows
		% Radius goes from 0 to 1 at edge, a little more in the corners.
		radius = sqrt((r - midY)^2 + (c - midX)^2) / min([midX, midY]);
		s(r, c) = min(1, radius); % Max out at 1
		h(r, c) = atan2d((r - midY), (c - midX));
	end
end

if threesixty == 0 
    
%%HSV wrap around
%     clear hsv
% h_circ = round(h*100+18000);
% cc = [hsv(18000); (hsv(18000))];

%equal colour space wrap around

% fullcolormatrix = importdata('colorwheel360.mat');
h_circ = round(h*2+360);

    cc = [fullcolormatrix(:,:); fullcolormatrix(:,:)]/255;

matric = zeros(270,270,3);

for col = 1:270
    for row = 1:270
        matric(col,row,:) = cc(h_circ(col,row),:);
    end
end
rgbImage = (matric); 
    
else
% Flip h right to left.
h = fliplr(mat2gray(h));
% Construct the hsv image.
jet = cat(3, h, s, v); %edit


% Construct the RGB image.
% rgbImage = hsv2rgb(jet);
% Display the RGB image.
% % % subplot(2, 2, 1:4);
% % % imshow(rgbImage, []);
% % % title('RGB Image, with V = 0.95', 'FontSize', fontSize);
% % % drawnow;

% Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
end
%% ring
% Create a logical image of a ring with specified
% inner diameter, outer diameter center, and image size.

% First create the image.
imageSizeX = size(rgbImage,1);%height-height*scaling;
imageSizeY = size(rgbImage,1);%height-height*scaling;
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = (height-height*scaling)/2;
centerY = (height-height*scaling)/2;
innerRadius = ((height-height*scaling)/2)*0.85;
outerRadius = (height-height*scaling)/2*0.98;
array2D = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2;
ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;
% ringPixels is a 2D "logical" array.
% Now, display it.
% % % % image(ringPixels) ;
% % % % colormap([0 0 0; 1 1 1]);
% % % % title('Binary Image of a Ring', 'Fo?hsvntSize', 25);



%% color circle

for i = 1:size(rgbImage,3)
    colourCircle(:,:,i) = rgbImage(:,:,i).* ringPixels*0.85;
    
%     colourCircle = matric.* ringPixels;
end

for x = 1:size(rgbImage,1)
  for y = 1:size(rgbImage,2)
    if ringPixels(x,y) == 0
      colourCircle(x,y,:) = 80/255;
    end
  end
end

cols = cc(1:2:720,:);
for i = 1:360
colourBar(i,:) = cols(mod(i+90,360)+1,:);
end