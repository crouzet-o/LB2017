% I've started creating a toolbox in order to be able to use functions
% that are defined outside a main script. This prevents one to have to
% define the same functions again each time a new script is generated.

% Aims of this script:
% - Extract formant frequencies for 2 paired CVC stimuli;
% - Generate a continuum of formant values / steps between the paired measurements;
% - Generate a continuum of formant transforms between the paired stimuli
% (2 versions: starting from each paired stimulus respectively);
% - Verify get_formants.m / get_formants.py 


addpath('~/work/lib/UMCG/LB2017/matlab/fmanipulation');

straight_path = '~/work/lib/UMCG/STRAIGHT/';
pyPraat_path  = '~/work/lib/UMCG/matlab/pyPraat/';

addpath(straight_path);






audiopath = 'recs-2017-03-27/Sem_Pip_01.flac';
[sound, fs] = audioread(audiopath);

% resampling to target sampling frequency
tfs = 44100;
tx = linspace(0, length(sound)/fs, length(sound));
sound = resample(sound, tx, tfs);
ofs = fs;
fs = tfs;

snd = audioplayer(sound, fs); 
%play(snd);


%t = (0:size(sp,2)-1)/1e3;

addpath(pyPraat_path);
% There's a weird issue with time output. Various options generate a
% different time vector
%[fx, r] = get_formants(audiopath, 'burg', 0.005, 5, 5500, 0.1);
fx = get_formants(audiopath);
rmpath(pyPraat_path);

% WARNING get_formants() operates on a file, whereas other aspects are
% processing the sound vector. And they don't have the same sampling
% frequency.

%plot(fx.t*1e3, fx.formants(:,1), '.');

% Remove all but the 3 first formants
keepF = 3;
fx.formants = fx.formants(:,1:keepF);
fx.bandwidths = fx.bandwidths(:,1:keepF);

% Interpolates formant frequencies over the original time vector (each
% sample is then associated with a measurement. Replaces all values outside
% the interval spanned by fx.t with NaNs (= NaN extrapolation, only 
% interpolation generates a value).
fx.formants = interp1(fx.t, fx.formants, tx, 'linear', NaN);


% fx.t contains time coordinates that correspond to the computed tracks,
% not to the sound wave. fx.time is the true time vector.
fx.time = tx;

% Interpolate F values
for i=1:size(fx.formants, 2)
    fx.formants(:,i) = interp_over_nan(fx.formants(:,i));
end

fx.bandwidths = fx.bandwidths(:,1:3);
fx.bandwidths = interp1(fx.t, fx.bandwidths, tx, 'linear', NaN);
for i=1:size(fx.bandwidths,2)
    fx.bandwidths(:,i) = interp_over_nan(fx.bandwidths(:,i));
end

smoothing = 1;
if smoothing~=0
    fx.formants   = smooth(fx.formants, smoothing);
    fx.bandwidths = smooth(fx.bandwidths, smoothing);
end



% plot(fx.time*1e3, fx.formants(:,1), '.');


%fx.formants(:, 1)
%fx.bandwidths(:, 1)
%fx.formants(floor(length(fx.formants(:, 1))/2), 1)
%fx.bandwidths(floor(length(fx.formants(:, 1))/2), 1)


[f0, ap] = exstraightsource(sound, fs);
sp = exstraightspec(sound, f0, fs);
f = (0:1024)/1024 * fs/2;
t = (0:length(f0)-1);
surf(t, f, 20 * log10(sp));
%set(h , 'zdata', get(h, 'xdata') * 0-1e3);

rgray = repmat([linspace(1, 0, 256)]', 1, 3);
colormap(rgray);

view(0,90);	
shading	flat;
%axis([0, max(t), 0, max(f)]);
xlim([0, max(t)]);
ylim([0, 5e3]);	
xlabel('Time (ms)');	
ylabel('Frequency (Hz)');


hold on;
%h = plot(tx*1e3, fx.formants, 'LineWidth', 3);	
% We use plot3() in order to control the z-axis and to prevent formant
% tracks to disappear under the spectrogram. A z-value of 1e6 is given to
% the formant points so it is always outperforming the spectrum energy.
h = plot3(tx*1e3, fx.formants, repmat(1e6, size(tx)), 'LineWidth', 3);	
hold off;



% STRAIGHT configuration
p = struct();

% VTL / f0 manipulation
p.Dvtl = 0;
p.Df0 = 0;

% Formant manipulation
p.smoothspec = 6; % number of spectral points (larger = smoother)
p.smoothing = 20; % Window-width in ms (integer)
p.DF1 = 3.5;
p.DF2 = -3.5;
p.DF3 = 0;

% manipulate_spectral_envelope.m applies on either a sound file or a
% corresponding vector. As it calls praat for formant extraction, it has to
% create a file from the vector.

% Verify necessity to estimate formant freq + use extended options (+pb maxNformants)?
% should we locate the vowel mid-position for formant measurement?
% if so, need for list of pairs + positions;
% anyway need list of pairs for generating continua;

addpath(pyPraat_path);
% There's a weird issue with time output. Various options generate a
% different time vector
%[fx, r] = get_formants(audiopath, 'burg', 0.005, 5, 5500, 0.1);
%fx = get_formants(audiopath);

[y, fs] = manipulate_spectral_envelope(sound, fs, p, true, true);
%[y, fs] = manipulate_spectral_envelope(sound, fs, p);
rmpath(pyPraat_path);

snd = audioplayer(sound, fs); 
play(snd);

pause(80/100)

snd = audioplayer(y, fs); 
play(snd);




% This example works fine. But I'm not sure it is okay for all pairs in the
% first Sem's recordings. I have to log information for all cases concerning original
% formant freqs (F1, F2, F3) and target freqs in order to identify possible
% issues. Also, it is based on changing F1 and F2 no matter what the value
% of the alternate in the pair would be (we just add / remove 3.5
% semitones).



