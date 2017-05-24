function [h, fx] = spectime(file, fparams)
  keepF = 5;
  audiopath = [file];
  [sound, fs] = audioread(audiopath);
  
  if fs > 44100
    tfs = 44100;
    tx = linspace(0, length(sound)/fs, length(sound));
    sound = resample(sound, tx, tfs);
  else
    tfs = fs;
    tx = linspace(0, length(sound)/fs, length(sound));
  end
  ofs = fs;
  fs = tfs;
  
  if nargin < 2
    fparams = struct();
  end
  
  if ~isfield(fparams, 'method')
    fparams.method = 'burg';
  end
  if ~isfield(fparams, 'time_step')
    fparams.time_step = 0.005;
  end
  if ~isfield(fparams, 'n_formants')
    fparams.n_formants = 5;
  end
  if ~isfield(fparams, 'max_freq')
      % max_freq should be validated based on sampling rate.
      % SR must be at least twice as large as max_freq.
      if fs < (2*5500)
        fparams.max_freq = floor(fs/2);
      else
        fparams.max_freq = 5500;
      end
  end
  if ~isfield(fparams, 'w_len')
      fparams.w_len = 0.1;
  end
  if ~isfield(fparams, 'export_method')
      fparams.export_method = 'matlabliteral';
  end
  
%  if fparams.w_len/fparams.time_step > length(sound)/fs
%    % w_len should must validated based on file duration. File needs to be
%    % at least twice as long as w_len.
%    % It also has to involve time_step. We need at least two time_steps
%    % within the sound duration.
%    fparams.w_len = floor(length(sound)/fs/2*100)/100
%    fparams.time_step = fparams.time_step/4
%  end

  
  % !!! One may probably favor to extract formants from a 
  % previously selected sound part in order to save processing time
  
  % There's a weird issue with time output. Various options generate a
  % different time vector
  %[fx, r] = get_formants(audiopath, 'burg', 0.005, 5, 5500, 0.1);
  [fx, r] = get_formants(audiopath, fparams.method, fparams.time_step, fparams.n_formants, fparams.max_freq, fparams.w_len, fparams.export_method);
  %fx = get_formants(audiopath);
  
  % WARNING get_formants() operates on a file, whereas other aspects are
  % processing the sound vector. And they don't have the same sampling
  % frequency.
  
  %plot(fx.t*1e3, fx.formants(:,1), '.');
  
  % Remove all but the 2 first formants
  fx.formants = fx.formants(:,1:keepF);
  %fx.bandwidths = fx.bandwidths(:,1:keepF);
  
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
  
  %fx.bandwidths = fx.bandwidths(:,1:3);
  %fx.bandwidths = interp1(fx.t, fx.bandwidths, tx, 'linear', NaN);
  %for i=1:size(fx.bandwidths,2)
  %    fx.bandwidths(:,i) = interp_over_nan(fx.bandwidths(:,i));
  %end
  
  smoothing = 1;
  if smoothing~=0
    fx.formants   = smooth(fx.formants, smoothing);
    %%fx.bandwidths = smooth(fx.bandwidths, smoothing);
  end
  
  
  
  [f0, ap] = exstraightsource(sound, fs);
  sp = exstraightspec(sound, f0, fs);
  f = (0:1024)/1024 * fs/2;
  t = (0:length(f0)-1);
  surf(t(1:size(sp, 2)), f, 20 * log10(sp));
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
  fname = strrep(file, '_', '\_');
  title(fname);

  
  hold on;
  %h = plot(tx*1e3, fx.formants, 'LineWidth', 3);	
  % We use plot3() in order to control the z-axis and to prevent formant
  % tracks to disappear under the spectrogram. A z-value of 1e6 is given to
  % the formant points so it is always outperforming the spectrum energy.
  h = plot3(tx*1e3, fx.formants, repmat(1e6, size(tx)), 'LineWidth', 3);
  hold off;
  
  % plot(fx.time*1e3, fx.formants(:,1), '.');
  
  
  %fx.formants(:, 1)
  %fx.bandwidths(:, 1)
end
