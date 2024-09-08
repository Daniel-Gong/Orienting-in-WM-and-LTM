%% Aedan Yue Li
% Memory & Perception Lab,
% University of Toronto
% September 25, 2018
%
% For questions, please email:
% aedanyue.li@mail.utoronto.ca
%% This script visualizes VCS space 
% REQUIRES: 
% 1. Psychtoolbox-3 (http://psychtoolbox.org/)
% 2. VCS space, download at: osf.io/d9gyf

%% This file should be placed in the same path as VCS stimuli: 
% define the path here -
PsychDefaultSetup(2);
path = 'VCS_space/VCSshapes/VCS_';

%% Define the sizing of the stimuli
Wheel.radius = 250; % Size of displayed wheel
stim.size = 300; % Size of presented stimulus

%% Presentation Logic
% Create window and setup params for where to show things:
Screen('Preference', 'SkipSyncTests', 1);
monitor = max(Screen('Screens'));

% Get the center of the monitor
[win, winRect] = Screen('OpenWindow', monitor);
  centerX = round(winRect(3)/2);
  centerY = round(winRect(4)/2);

  % define the size of presented wheel and stimuli based on defined values
  Wheel.rect = CenterRect([0 0 Wheel.radius*2 Wheel.radius*2], winRect);
  stimRect = CenterRect([0 0 stim.size stim.size], winRect);
  Screen('Flip', win);

% Logic for the shape wheel
mouse_click = [];
oldAngle = -1;
    while ~any(mouse_click)
      [curX,curY, ~] = GetMouse(win); % you can end the while loop by 
      % replacing '~' with "mouse_click", so that clicking a shape can 
      % start the next trial. However, you would need to place this code 
      % within another loop
      curAngle = GetPolarCoordinates(curX,curY,centerX,centerY); % get the 
      % current angle the GetPolarCoordinates function, taking the current 
      % mouse xy position, and the centre of the screen
      
      % this part draws the position of the mouse cursor on the wheel,
      % using the polar2xy function
      [dotX1, dotY1] = polar2xy(curAngle,Wheel.radius-5,centerX,centerY);
      [dotX2, dotY2] = polar2xy(curAngle,Wheel.radius+5,centerX,centerY);
      
      % Draw everything
      Screen('FrameOval', win, [128,128,128], Wheel.rect);
      Screen('DrawLine', win, [0 0 0], dotX1, dotY1, dotX2, dotY2, 4);
      
      %% Now display the shapes
      % If angle changed, close old texture and make new one in correct color:
      if (curAngle ~= oldAngle) && round(curAngle) ~= 0 % at the start of the trial, oldAngle == -1 from previous code
          
        % calculate angular position from mouse x y coordinates, then load
        % the associated texture from stimuli path
        curAngRounded = round(curAngle);
        shape_angle = curAngRounded;
        newImg = imread([path num2str(shape_angle) '.jpg']);
        newImg = imresize(newImg, [stim.size stim.size]);
        curTexture = Screen('MakeTexture', win, newImg);
      end

      % Show stimulus:
      Screen('DrawTexture', win, curTexture, [], stimRect);
      Screen('Flip', win);
      oldAngle = curAngle;
      
      % Allow user to quit on each frame:
      [~,~,keys]=KbCheck;
      if keys(KbName('q'))
        sca; error('User quit');
      end
    end
        
% Determine  angle depending on the position of cursor relative to the
% center of the screen   
function [angle, radius] = GetPolarCoordinates(h,v,centerH,centerV)
  % get polar coordinates
  hdist   = h-centerH;
  vdist   = v-centerV;
  radius     = sqrt(hdist.*hdist + vdist.*vdist)+eps;
  
  % determine angle using cosine (hyp will never be zero)
  angle = acos(hdist./radius)./pi*180;
  
  % correct angle depending on quadrant
  angle(hdist == 0 & vdist > 0) = 90;
  angle(hdist == 0 & vdist < 0) = 270;
  angle(vdist == 0 & hdist > 0) = 0;
  angle(vdist == 0 & hdist < 0) = 180;
  angle(hdist < 0 & vdist < 0)=360-angle(hdist < 0 & vdist < 0);
  angle(hdist > 0 & vdist < 0)=360-angle(hdist > 0 & vdist < 0);
end

% given the current angle of the mouse cursor, radius of the circle and
% the center of the screen, define the x and y coordinates  
function [x, y] = polar2xy(angle,radius,centerH,centerV)  
  x = round(centerH + radius.*cosd(angle));
  y = round(centerV + radius.*sind(angle));
end