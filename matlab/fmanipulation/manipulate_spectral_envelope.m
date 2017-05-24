function [y, fs] = manipulate_spectral_envelope(varargin)
% TODO: write a documentation


% Check if input is a filename or sound vector
if ischar(varargin{1})
    wavIn = varargin{1};
    i_arg = 2;
    inf = audioinfo(wavIn);
    if inf.NumChannels~=1
        error('Can only process single channel audio (%d channel(s) found). To process multichannel audio files, separate the channels to individual files or vectors.', inf.NumChannels);
    end
elseif isnumeric(varargin{1})
    if ~isnumeric(varargin{2})
        error('If the first argument is a vector, the second argument must be the sampling frequency');
    end
    x  = varargin{1};
    fs = varargin{2};
    i_arg = 3;
    if size(x,2)~=1
        error('Can only process single channel audio (%d channel(s) found). To process multichannel signals, separate the channels to individual vectors.', size(x,2));
    end
    wavIn = [tempname(), '.wav'];
    audiowrite(wavIn, x, fs);
end

% Parameters
p = varargin{i_arg};

% Force rebuild everything
if length(varargin)<i_arg+1
    force_rebuild = false;
else
    force_rebuild = varargin{i_arg+1};
end

% Debug
if length(varargin)<i_arg+2
    debug = false;
else
    debug = varargin{i_arg+2};
end


%----------
% Do we already have the sound we want?

[d, b]  = fileparts(wavIn);
matFile = fullfile(d, [b, '.straight.mat']);
wavOut  = [make_fname(wavIn, p), '.wav'];

if exist(wavOut, 'file') && ~force_rebuild
    if debug
        fprintf('A file has been found that seems to fit the requested parameters: "%s". It will be loaded.\n', wavOut);
    end
    [y, fs] = audioread(wavOut);
    return
end

%----------
% STRAIGHT analysis

if debug
    fprintf('STRAIGHT analysis for "%s":\n', wavIn);
end
if exist(matFile, 'file')
    if debug
        fprintf('  Loading mat file "%s".\n', matFile);
    end
    load(matFile);
else
    [x, fs] = audioread(wavIn);

    if debug
        fprintf('  Extracting source.\n');
    end
    [f0, ap] = exstraightsource(x, fs);

    if debug
        fprintf('  Extracting spectral envelope.\n');
    end
    sp = exstraightspec(x, f0, fs);
    x_rms = rms(x);

    if debug
        fprintf('  Saving to mat file "%s".\n', matFile);
    end
    save(matFile, 'fs', 'f0', 'sp', 'ap', 'x_rms');
end

if size(sp,1)<1025
    error('The sampling frequency of the input file (%d Hz) is insufficient.', fs);
end

t = (0:size(sp,2)-1)/1e3;

%----------
% Manipulation parameters

% SER and VTL
if isfield(p, 'Dvtl')
    ser = 2^(-p.Dvtl/12);
elseif isfield(p, 'ser')
    ser = p.ser;
else
    ser = 1.0;
end

% F0
if isfield(p, 'f0')
    f0(f0~=0) = p.f0;
elseif isfield(p, 'mf0')
    f0 = f0/exp(mean(log(f0(f0~=0))))*p.mf0;
elseif isfield(p, 'Df0')
    f0 = f0*2^(p.Df0/12);
end

% Formants
DF = ones(1,3);
calculate_formants = true;
for i=1:3
    Df = sprintf('DF%d', i);
    Hf = sprintf('HF%d', i);
    if isfield(p, Df)
        DF(i) = 2^(p.(Df)/12);
        DF_type{i} = 'relative'; % As a ratio
        if DF(i)~=1
            calculate_formants = true;
        else
            DF_type{i} = 'none';
        end
    elseif isfield(p, Hf)
        DF(i) = p.(Hf);
        DF_type{i} = 'absolute'; % In Hertz
        if DF(i)~=0
            calculate_formants = true;
        else
            DF_type{i} = 'none';
        end
    else
        DF_type{i} = 'none';
    end
end

if debug
    disp(DF);
end

%Smoothing
if isfield(p, 'smoothing')
    % Smoothing of the formant values, window width in ms (integer)
    smoothing = p.smoothing;
else
    smoothing = 0;
end

% Added by OC to smooth spectral response
if isfield(p, 'smoothspec')
    % Smoothing of the spectral envelope, window width (integer) in number 
    % of spectral points
    if p.smoothspec < 3 && p.smoothspec > 0
        smoothspec = 3;
    else
        smoothspec = round(p.smoothspec);
    end
else
    smoothspec = 0;
end

%----------
% Get formants if necessary

if calculate_formants
    if debug
        fprintf('  Calculating formants...\n');
    end
    formants = get_formants(wavIn);
    
    % Resample the formant matrix to the right time scale and retaining
    % only the first three formants
    formants.formants = formants.formants(:,1:3);
    formants.formants = interp1(formants.t, formants.formants, t, 'linear', NaN);
    for i=1:size(formants.formants,2)
        formants.formants(:,i) = interp_over_nan(formants.formants(:,i));
    end
    
    formants.bandwidths = formants.bandwidths(:,1:3);
    formants.bandwidths = interp1(formants.t, formants.bandwidths, t, 'linear', NaN);
    for i=1:size(formants.bandwidths,2)
        formants.bandwidths(:,i) = interp_over_nan(formants.bandwidths(:,i));
    end
    
    if smoothing~=0
        formants.formants   = smooth(formants.formants, smoothing);
        formants.bandwidths = smooth(formants.bandwidths, smoothing);
    end

    %----------
    % Spectral envelope manipulation (if necessary)
    
    % NOTE: we are not applying the formant changes to AP, but maybe we
    % should...
    
    if debug
        fprintf('  Manipulating spectral envelope based on formant frequencies...\n');
    end
    f = (0:1024)/1024*fs/2;
    
    sp_new = zeros(size(sp));
    
    for it=1:length(t)
        anchor_ref = 0;
        anchor_new = 0;
        for i=1:3
            fp = [-1, 0, 1] * formants.bandwidths(it, i) / 2 + formants.formants(it, i); % Here we are assuming that there will be not overlap between formants +/- bandwith...
            anchor_ref = [anchor_ref, fp];
            switch DF_type{i}
                case 'relative'
                    fp_new = DF(i) .^ [1, 1, 1] .* fp;
                case 'absolute'
                    fp_new = DF(i) * [1, 1, 1] + fp;
                case 'none'
                    fp_new = fp;
            end
            anchor_new = [anchor_new, fp_new];
        end
        anchor_ref = [anchor_ref, f(end)];
        anchor_new = [anchor_new, f(end)];
        
        if debug && it==350
            figure(203)
            plot([0, fs/2], [0, fs/2], '--', 'color', [1 1 1]*.7)
            hold on
            plot(anchor_ref, anchor_new, '-+')
            plot(anchor_ref(3:3:9), anchor_new(3:3:9), 'o')
            hold off
            xlim([0, 5e3]), ylim([0, 5e3]);
            xlabel('Input frequency (Hz)');
            ylabel('Output frequency (Hz)');
        end
        
        c = interp1(anchor_new, anchor_ref, f, 'linear');
        
        sp_new(:,it) = interp1(f, sp(:,it), c);

        %size(sp_new)
        %if debug && it==350
        %    figure(101)
        %    hold('on'); plot(sp_new); hold('off');
        %end
        %if debug && it==350
        %    figure(101);
        %    plot(sp_new);
        %end
        %size(sp_new(:, it))
        if smoothspec ~= 0
            w = hann(smoothspec);
            w = w / sum(w);
            xx = sp_new(:, it);
            xx = conv2([ones(smoothspec,1) * xx(1); xx(:); ones(smoothspec,1) * xx(end)], w, 'same');
            xx = xx(smoothspec+1:end-smoothspec);
            sp_new(:, it) = xx;
        end
         
        if debug && it==350
            figure(223)
            h = plot(f, sp(:,it));
            hold on
            h_new = plot(f, sp_new(:, it));
            plot(formants.formants(it, :), interp1(f, sp(:,it), formants.formants(it, :)), 'v', 'Color', h.Color);
            plot(anchor_new(3:3:9), interp1(f, sp_new(:,it), anchor_new(3:3:9)), '+', 'Color', h_new.Color);
            hold off
            xlim([0, 5000]), set(gca, 'YScale', 'log'), xlabel('Frequency (Hz)');
        end
    end
    
    if debug
        figure(198)
        subplot(1,2,1)
        pcolor(sp)
        shading flat
        subplot(1,2,2)
        pcolor(sp_new)
        shading flat
        
        figure(203)
        hold off
    end
else
    sp_new = sp;
end

straight_p.frequencyAxisMappingTable = ser;

if debug
    fprintf('  Re-synthesizing...\n');
end
y = exstraightsynth(f0, sp_new, ap, fs, straight_p);

y = y/rms(y)*x_rms;
if max(abs(y))>1
    warning('Output was renormalized for "%s" (peak at %f).', wavOut, max(abs(y)));
    y = 0.98*y/max(abs(y));
end

audiowrite(wavOut, y, fs);

