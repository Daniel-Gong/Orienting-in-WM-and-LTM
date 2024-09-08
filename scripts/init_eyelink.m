function  init_eyelink(screenX, screenY)

%%=== CREDITS ===
%Author: Gianpiero Monittola @ University of Trento
%Date: June 2019
%Last Version: 11.09.2019
%MATLAB R2012b
%Psychtoolbox-3-Psychtoolbox-3.0.11_for_MATLAB32bits

% === Script Description ===
%Force default setting when the host PC (re)boots.

     %Set Camera.Tracking mode
     Eyelink('Command', 'sample_rate = 1000'); % Select the sampling rate for the recording. Here 1000 Hz is select
     Eyelink('Command', 'use_ellipse_fitter = NO'); % Select the tracking mode (Pupil-only vs. Pupil-CR)
     
     Eyelink('Command', 'enable_search_limits = ON'); % Autothreshold on mouse click on setp mode image
     Eyelink('Command', 'autothreshold_click = TRUE');
     Eyelink('command', 'elcl_tt_power = %d',1); % set illumination power in camera setup screen 1 = 100%, 2 = 75%, 3 = 50%
     
     %Set Option.Calibration and Validation
     Eyelink('command', 'calibration_type = HV5'); % HV9 = 9 points, biquadratic with corner correction H3. HV3. HV5. HV9. HV13
     Eyelink('Command', 'generate_default_targets = NO');
     Eyelink('Command', 'automatic_calibration_pacing = 1000');	% delay in ms between calibratio and valedation
	 Eyelink('Command', 'randomize_calibration_order = YES'); % randomize the calibration and validation target presentation order
     Eyelink('Command', 'cal_repeat_first_target = YES'); % redisplay the fist calibration and validation target at the end of the calibration sequence 
	 Eyelink('Command', 'val_repeat_first_target = YES');
     
     %Set Option.Tracking
     Eyelink('Command', 'enable_search_limits = YES'); %enable display of global search limit
     Eyelink('Command', 'track_search_limits = NO'); %enable tracking of pupil to global search limit
     Eyelink('command', 'pupil_size_diameter = AREA'); %select the type of data used for pupil size
     
     %Set Option.Event and Data Processing
     Eyelink('Command', 'recording_parse_type = GAZE'); %data type used to compute velocity for parsing of eye movement during recording
     Eyelink('Command', 'select_parser_configuration = 0'); %select the preset standard parser setup
     Eyelink('Command', 'heuristic_filter = ON'); 
     Eyelink('Command', 'heuristic_filter = 1 2'); %1 File Sample Filter EXTRA, 2 Link/Analog Filter STD
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % 1 File Sample Filter EXTRA
     % EyeLink use a heuristic filtering algorithm for data smoothing
     % Data filtering can be applied indipendently for the data saved in
     % the EDF file e for data set to link/analog output.
     % The current option selects filter level of data recorded to the EDF file.
     % Each increase in filter level reduces noise by factor of 2 to 3.
     
     % 2 Link/Analog Filter STD
     % Select the filter level for data available via Ethernet link and analog output 
     % Each increase in filter level reduces noise by factor of 2 to 3 
     % but introduces a 1-sample delay to the link sample feed.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % from chage demo (psychtoolbox) ####start
     % set parser (conservative saccade thresholds)
     Eyelink('command', 'saccade_velocity_threshold = 35');
     Eyelink('command', 'saccade_acceleration_threshold = 9500');
     
     % This command is crucial to map the gaze positions from the tracker to
     % screen pixel positions to determine fixation
     Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenX-1, screenY-1);
     
     Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenX-1, screenY-1);
     % ####end
     
     
     %% Set Option.Analog Output
     % Select the type of the analog output  
     Eyelink('Command', 'analog_out_data_type = GAZE'); % GAZE is scren gaze x,y 
     % RAW is uncalibrated pupil x,y in camera coordinate;
     % HREF is head referenced calibrated x,y
     
      
     %Set Option.File data Contents
     Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT,FIXUPDATE');
     Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
     
     % This commad sets the contents of the sample data in the EDF file recording
     % Selecting ‘Samples’ will record data samples to the EDF file, and selecting Events will record online parsed events.
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % The data contents of an EDF file are organized in two streams: samples and events. 
     % Samples are used to record instantaneous eye position data, while events are used to record important occurrences, 
     % either from the experimental application or from changes in the eye data.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Set Option.File Sample Contents
     Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT,FIXUPDATE');
     Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % The data passed are in a list with these types:
     % LEFT,RIGHT: data for one or both eyes
     % GAZE: screen xy (gaze) position (pupli position for calibration)
     % GAZERES: units-per-degree screen resolution (start end of the event)
     % HREF: head-referenced gaze position
     % PUPIL: raw eye camera pupil coordinates
     % AREA: pupil area
     % VELOCITY: velocity of parsed position-type (avg, peak)
     % Status: warning and error flags, aggregated accros event
     % FIXAVG: includE only averages in fixation and events, to reduce file size
     % NOSTART: start events have no data, just time stamp
     % BUTTON: button state and change flags
     % INPUT: input port data lines
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Set Option.Recording data view
     Eyelink('Command', 'rec_plot_enabled = NO');
     Eyelink('Command', 'rec_plot_data = GAZE'); % controls what to show on the Record screen during data output. 

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % If Record View is set to Gaze Cursor, the Host PC Record screen will display the participant’s current gaze
     % position as a cursor graphic overlaid on a simulated display screen. If set to Plotting, x, y data traces 
     % will be graphed as a function of time. The user can further select which data type should be plotted.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     %Set Option.Video overlay  
     Eyelink('Command', 'video_overlay_on = OFF'); %sets the state of the video overlay mode. 
     %This is updated at the end of each session to the lastrun.ini file,
     %which overrides the vidvl.ini file. OFF if overy mode is off.
     %Clicking ‘Video Setup’ goes to the Video Setup screen. Clicking ‘Enable Overlay’ activates the video overlay option.

     % calibration/drift correction target
     Eyelink('Command', 'button_function 5 "accept_target_fixation"');
     
end