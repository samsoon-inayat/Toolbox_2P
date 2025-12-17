function animal = process_behavior_signals(animal)

if ~exist('animal','var')
    animal = evalin('base','animal')
end
for an = 1:numel(animal)
    % filename =
    tanim = animal(an);
    videos = tanim.video;
    v = VideoReader("your_video.mp4");
    fps = v.FrameRate;   % should be ~60
    
    % Read first frame and pick a circular ROI
    frame1 = readFrame(v);
    imshow(frame1);
    title("Draw circular ROI over LED, double-click to finish");
    h = drawcircle;      % requires newer MATLAB versions
    mask = createMask(h);
    
    % Reset reader
    v.CurrentTime = 0;
    
    led_mean = [];
    while hasFrame(v)
        fr = readFrame(v);
        g = rgb2gray(fr);
        led_mean(end+1,1) = mean(g(mask));
    end
    
    t = (0:length(led_mean)-1)'/fps;
    
    % Auto threshold (percentile midpoint)
    lo = prctile(led_mean, 10);
    hi = prctile(led_mean, 90);
    thr = (lo+hi)/2;
    
    air_on = led_mean >= thr;
    
    figure; plot(t, led_mean); yline(thr,'--'); title("LED ROI mean intensity");
    figure; stairs(t, air_on); ylim([-0.2 1.2]); title("Air ON (0/1)");
    save("air_pulse_from_led.mat","t","led_mean","air_on","fps","thr");
end