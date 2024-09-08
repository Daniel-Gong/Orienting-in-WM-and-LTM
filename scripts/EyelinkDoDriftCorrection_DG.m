function [success,repeat] = EyelinkDoDriftCorrection_DG(el, x, y, draw, allowsetup)
%% This is a modified version of the original drift correction function
%%
success=0;
% DO PRE-TRIAL DRIFT CORRECTION */
% We repeat if ESC key pressed to do setup. */
% Setup might also have erased any pre-drawn graphics. */

% if no x and y are supplied, set x,y to center coordinates
if ~exist('x', 'var') || isempty(x) || ~exist('y', 'var') || isempty(y)
	[x,y] = WindowCenter(el.window); % convenience routine part of eyelink toolbox
end

if ~exist('draw', 'var') || isempty(draw)
	draw=1;
end

if ~exist('allowsetup', 'var') || isempty(allowsetup)
	allowsetup=1;
end

repeat = 0;
while 1
    if Eyelink('IsConnected')~=1   % Check link often so we don't lock up if tracker lost
        success=-1;
        sca;
        return;
    end
    % DRIFT CORRECTION */
    % 3rd argument would be 0 to NOT draw a target */
    % fprintf('drifcorr at % d %d\n', x, y );
    error = EyelinkDoDriftCorrect_DG(el, x, y, draw, allowsetup);
    HideCursor;
    if error==el.TERMINATE_KEY
        success=0;
        return;
    end
    % repeat if ESC was pressed to access Setup menu
    if error==el.ESC_KEY && repeat <= 1
        repeat = repeat + 1;
    elseif error==el.ESC_KEY && repeat > 1
        success=-1;
        return; 
    end        
    if error~=el.ESC_KEY
        success=1;
        return; 
    end
end
end
