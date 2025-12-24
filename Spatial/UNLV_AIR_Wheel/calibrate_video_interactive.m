function calibrate_video_interactive(video_path)
   

vp = evalin('base','vp');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');
n = 0;
%%
% 2. Read the first frame
firstFrame = readFrame(vp);

% 3. Display the frame
figure('Name', 'Calibration Tool', 'Position', [100, 100, size(firstFrame, 2)/2, size(firstFrame, 1)/2]);
imshow(firstFrame);
title({'Step 1: Draw a line along a known distance.'; 'Step 2: Press Enter in the Command Window.'});

% 4. Use imdistline to measure a distance interactively
h = imdistline(gca); 

% 5. Wait for user to confirm in the command window
fprintf('\n***ACTION REQUIRED: Draw your line on the figure, then press the Enter key in the command window when ready.***\n');
pause; % Wait for user to press any key/Enter

% 6. Get the distance once the user is ready
api = iptgetapi(h);
pixelDistance = api.getDistance();

% 7. Prompt the user for the actual distance and units
prompt = {'Enter the known distance:', 'Enter the units (e.g., mm, cm):'};
dlgtitle = 'Real-World Calibration Input';
dims = [1 35];
definput = {'10', 'cm'}; % Default values
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('Calibration cancelled.');
    scaleFactor = [];
    return;
end

knownDistance = str2double(answer{1});
units = answer{2};

% 8. Calculate the scale factor (e.g., 0.073 cm/pixel)
scaleFactor = knownDistance / pixelDistance;

% 9. Display the final results in the command window
fprintf('\n--- Calibration Results ---\n');
fprintf('Pixel distance: %.2f pixels\n', pixelDistance);
fprintf('Known distance: %.2f %s\n', knownDistance, units);
fprintf('Scale Factor: %.4f %s/pixel\n', scaleFactor, units);
fprintf('---------------------------\n');

% Close the figure window automatically
close(gcf);