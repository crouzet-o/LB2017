tic;

nbpoints = 7;
keepF = 3;

debug = false;
wavestorage = true;
if debug
    wavestorage = true;
end
if wavestorage
    savedir = 'base-sounds/';
    status = mkdir(savedir);
end

logfile = 'log-soundprocessing-basesounds.csv';
%csvsep = ';';
logid = fopen(logfile, 'w');
fprintf(logid, '%s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s ; %s\n', 'output', 'startpoint', 'endpoint', 'numContinuum', 'startF1', 'endF1', 'targetF1', 'zF1', 'stF1', 'startF2', 'endF2', 'targetF2', 'zF2', 'stF2');


addpath('~/work/lib/UMCG/matlab/fmanipulation');

straight_path = '~/work/lib/UMCG/STRAIGHT/';
pyPraat_path  = '~/work/lib/UMCG/matlab/pyPraat/';
addpath(straight_path);
addpath(pyPraat_path);

% Load CVC audio files pairings

%fid = fopen('pairings-full.csv');
%C = textscan(fid, '%n %s %s %s %s %s %s %s %n %n %n %n', 'Delimiter', ',', 'TreatAsEmpty', {'NA', 'na', '--undefined--'}, 'EmptyValue', NaN);
%fclose(fid);

% celldisp(C)

%contrast = C{2};
%axis = C{3};
%cXc = C{4};

% Load time positions for middle of V_1 and V_2 (TODO)
%ex1mid = C{11};
%ex2mid = C{12};

%ex1 = strcat(C{5},'_0',num2str(C{9}),'.flac');
%ex2 = strcat(C{6},'_0', num2str(C{10}),'.flac');
%clear('contrast','axis','cXc', 'ex1mid', 'ex2mid', 'ex1', 'ex2');


% Declare speakers' names (soundfile prefixes)
speakers = {'Sem'};

% Load main stimulus file
fid = fopen('pairings.csv');
%fid = fopen('pairings_test.csv');
C = textscan(fid, '%s %s %s %s %s %s %s %n %n', 'Delimiter', ';', 'HeaderLines', 1, 'TreatAsEmpty', {'NA', 'na', '--undefined--'}, 'EmptyValue', NaN);
fclose(fid);

nstim = length(C{1});

% Generate random file pairings
%rng('shuffle');
% Fix the random seed in order to get the same results each time the random
% permutation is performed
%rng(pi);

pairs1 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
pairs2 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
contrast = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
axis = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
cXc = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
w1 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
w2 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
ex1 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);
ex2 = mat2cell(repmat(NaN, length(C{1}), 1), length(C{1}), 1);

npairs = 4;
for k = 1:nstim
    rng(k*pi); % Set the random seed to select the first member of the pair
    pairs1{k} = randperm(npairs)';
    rng(length(C{1})+(k*pi)) % Set the random seed to select the second member of the pair
    pairs2{k} = randperm(npairs)';

    contrast{k} = repmat(C{1}(k), npairs, 1);
    axis{k} = repmat(C{2}(k), npairs, 1);
    cXc{k} = repmat(C{3}(k), npairs, 1);
    
    w1{k} = repmat(C{4}(k), 4, 1);
    w2{k} = repmat(C{5}(k), 4, 1);
    ex1{k} = strcat(repmat(C{4}(k),4,1), '_0', num2str(pairs1{k}), '.flac');
    ex2{k} = strcat(repmat(C{5}(k),4,1), '_0', num2str(pairs2{k}), '.flac');

end

% Load time points for formant extraction from companion file
fid = fopen('midV.csv');
tP = textscan(fid, '%s %s %n %n %n %n', 'Delimiter', ';', 'HeaderLines', 1, 'TreatAsEmpty', {'NA', 'na', '--undefined--'}, 'EmptyValue', NaN);
fclose(fid);

timepoints = cell2mat(tP(3:(3+npairs-1)));

%%

for s = 1:length(speakers) % For each speaker
    for i = 1:nstim % For each contrast available
        for j = 1:npairs % For each exemplar of a pair
            stim = {cell2mat(strcat(speakers{s}, '_', ex1{i}(j))), cell2mat(strcat(speakers{s}, '_', ex2{i}(j)))};
            
            ftargets = nan(2, keepF);
            for k = 1:2 % /start -> end/ then /end -> start/
                % Two versions are generated depending on whether we start from
                % one or the alternate endpoint
                
                % Compute formant freqs for the endpoints
                
                % Define stimstart / stimend + initialize storage for formant endpoints
                kcomp = setxor(k, 1:2);
                stimstart = stim{k};
                stimend = stim{kcomp};
                stimseq = {stimstart, stimend}
                %if (debug == 1)
                %    ['Start: ', stimstart, ' End: ', stimend]
                %end
                
                
                % Determination of file / midpoint correspondence
                
                % The strfind() call is performed on a string, therefore it
                % returns a string. It has to be converted to numeric then as
                % we use it to generate array-indexing.
                exN = cellfun(@str2num, {stimseq{1}(strfind(stimseq{1}, '.flac')-1), stimseq{2}(strfind(stimseq{2}, '.flac')-1)});
                
                midpoints = cell2mat({timepoints(find(strcmp(tP{1}, speakers{s}) & strcmp(tP{2}, w1{i}{1})), exN(1)),...
                    timepoints(find(strcmp(tP{1}, speakers{s}) & strcmp(tP{2}, w2{i}{2})), exN(2))});
                
%            end            
% Alternative approach
% Take the pair
% Extract formants first for each endpoint
% Store information
% Loop through each pair and process sound

% We've ended the loop above in order to store the relevant double-entry
% data in a array (one entry for startstim, one for endstim) and we will 
% simply reverse / fliplr() the array when going reverse
%                for k = 1:2
                    % Extract formants for start (l=1) and end (l=2) points
                    
                    % Formant extraction at required time points
                    audiopath = ['recs-2017-03-27/', stimstart];
                    [sound, fs] = audioread(audiopath);
                    
                    tfs = 44100;
                    tx = linspace(0, length(sound)/fs, length(sound));
                    if (fs ~= tfs)
                        sound = resample(sound, tx, tfs);
                        ofs = fs;
                        fs = tfs;
                    end
                    % !!! One may probably favor to extract formants from a 
                    % previously selected sound part in order to save processing time
                    
                    % Better generate it depending on speaker !!!!
                    fparams = struct('method', 'burg', 'time_step', 0.005, 'n_formants', 6, 'max_freq', 5500, 'w_len', 0.05, 'export_method', 'matlabliteral');
                    fx = get_formants(audiopath, fparams.method, fparams.time_step, fparams.n_formants, fparams.max_freq, fparams.w_len, fparams.export_method);
                    
                    % Remove all but the 2/3 first formants
                    fx.formants = fx.formants(:, 1:keepF);
                    %fx.bandwidths = fx.bandwidths(:, 1:keepF);
                    
                    % Interpolates formant frequencies over the original time vector (each
                    % sample is then associated with a measurement. Replaces all values outside
                    % the interval spanned by fx.t with NaNs (= NaN extrapolation, only 
                    % interpolation generates a value).
                    fx.formants = interp1(fx.t, fx.formants, tx, 'linear', NaN);
                    
                    % fx.t contains time coordinates that correspond to the computed tracks,
                    % not to the sound wave. fx.time is the true time vector.
                    fx.time = tx;
                    
                    % Interpolate F values
                    for f = 1:size(fx.formants, 2)
                        fx.formants(:, f) = interp_over_nan(fx.formants(:, f));
                    end
                    
                    %fx.bandwidths = fx.bandwidths(:,1:3);
                    %fx.bandwidths = interp1(fx.t, fx.bandwidths, tx, 'linear', NaN);
                    %for f = 1:size(fx.bandwidths,2)
                    %    fx.bandwidths(:, f) = interp_over_nan(fx.bandwidths(:, f));
                    %end
                    
                    smoothing = 1;
                    if smoothing ~= 0
                        fx.formants   = smooth(fx.formants, smoothing);
                        %fx.bandwidths = smooth(fx.bandwidths, smoothing);
                    end
                    
                    % Use time position -> convert to index
                    %pos = round(midpoints(l)/1000*fs/1); % measured position
                    % 30ms around time position -> converted to index
                    window = 30;
                    pos = [round((midpoints(k)-(window/2))/1000*fs/1):round((midpoints(k)+(window/2))/1000*fs/1)];
                    %pos = floor(2 * length(fx.formants(:, 1))/4) % V-mid
                    freqstraj = fx.formants(pos, :);
                    if length(pos)>1
                        freqs = mean(freqstraj, 1); % Compute the mean freq over the array lines (1st dimension)
                    end
                    %trajcenter = fx.formants((floor(1 * length(fx.formants(:, 1))/4)):(floor(3 * length(fx.formants(:, 1))/4)), :); % Middle-zone
                    %mFormant = median(trajcenter, 1);
                    ftargets(k,:) = freqs;

            end
            % Generate a linear continuum of formant values between V_1 and
            % V_2 (and between V_2 and V_1)
            
            % First compute the continuum
            ztargets = f2Bark(ftargets);
            
            % We define the euclidean distance vector between the two 
            % endpoints and generate the vector of required euclidean 
            % distances. Note that euclidist() will compute the Euclidean 
            % distance between 2 points for any number of dimensions (2, 
            % 3, 4...).
            ed = euclidist(ztargets(1,:), ztargets(2,:)); % fmanipulation toolbox
            dvect = linspace(0, ed, nbpoints);
            
            % Then we compute the y = ax + b equivalent function between the 
            % two endpoints. This is where introducing a 3rd dimension
            % si challenging...
            [sl, yi] = affine(ztargets(2,:), ztargets(1,:)); % fmanipulation toolbox
                
                % We generate the sequence of linearly spaced
                % points between the endpoints
                if (l == 1)
                    target_dvect = dvect;
                    target_sl = sl;
                else
                    target_dvect = dvect;
                    target_sl = sl;
                end
                
                [coords, theta] = lpoints(ztargets(1,:), target_dvect, target_sl, [1,2]); % fmanipulation toolbox
                %   dvect, -sl good for F1, but -dvect, sl good for F2
                    
                % A simpler and equivalent way to do the continuum computation 
                % is to request a linear sequence of F1, F2, F3 values
                % independently. This is equivalent to a Euclidean space
                % computation because we generate the continuum on a linear
                % equation
                if (keepF == 2)
                    coords = [linspace(ztargets(1,1), ztargets(2,1), nbpoints) ;...
                        linspace(ztargets(1, 2), ztargets(2, 2), nbpoints)];
                elseif (keepF == 3)
                    coords = [linspace(ztargets(1,1), ztargets(2,1), nbpoints) ;...
                        linspace(ztargets(1, 2), ztargets(2, 2), nbpoints) ;...
                        linspace(ztargets(1, 3), ztargets(2, 3), nbpoints)];
                else % In any other case we only keep 1 formant
                    coords = linspace(ztargets(1,1), ztargets(2,1), nbpoints);
                end

                % The Bark coordinates are then converted back to Hertz coordinates
                tfreqs = f2Bark(coords, 'back');

                for l = 1:2
                    % Generate the continuum (2 versions starting from the 2 endpoints)
                                    
                    % Define stimstart / stimend + initialize storage for formant endpoints
                    lcomp = setxor(l, 1:2);
                    stimstart = stim{l};
                    stimend = stim{lcomp};
                    stimseq = {stimstart, stimend}
                    if (debug == 1)
                        ['Start: ', stimstart, ' End: ', stimend]
                    end

                    if (l == 2)
                        tfreqs = fliplr(tfreqs);
                        coords = fliplr(coords);
                    end
                    
                    %if (l == 1)
                    %    %if (ztargets(1(l) > ztargets{1}(setxor(l, 1:2)))
                    %    target_dvect = dvect;
                    %else
                    %    target_dvect = -dvect;
                    %end
                    %    target_dvect = -dvect;
                    
                    %if  (ztargets{2}(l) > ztargets{2}(setxor(l, 1:2)))
                    %target_sl = sl;
                    %else
                    %    target_sl = -sl;
                    %end
                    % It should probably test positions within the formant space... not
                    % order of processing because we may start with I-i for
                    % example in some pairs...
                    % For the A->B continuum
                            
                    audiopath = ['recs-2017-03-27/', stimseq{l}];
                    [sound, fs] = audioread(audiopath);
                    
                    tfs = 44100;
                    tx = linspace(0, length(sound)/fs, length(sound));
                    sound = resample(sound, tx, tfs);
                    ofs = fs;
                    fs = tfs;
                    
                    [f0, ap] = exstraightsource(sound, fs);
                    sp = exstraightspec(sound, f0, fs);
                    
                    % STRAIGHT configuration
                    p = struct();
                    
                    % VTL / f0 manipulation
                    p.Dvtl = 0;
                    p.Df0 = 0;
                    
                    % Formant manipulation
                    % !!! Convert from tfreqs in Barks to semitones relation
                    % with ref = tfreqs(1) or tfreqs(length(tfreqs))
                    p.smoothspec = 6; % number of spectral points (larger = smoother)
                    p.smoothing = 20; % Window-width in ms (integer)
                    
                    %if (l == 1)
                    %    tfreqs = tfreqs;
                    %else
                    %    tfreqs = fliplr(tfreqs);
                    %end
                    
                    for m = 1:(size(tfreqs, 2))
                        %ftargets{1}
                        %ftargets{2}
                        %tfreqs
                        output = [stimstart, '_to_', stimend, '_', num2str(m-1), '.wav'];
                        p.DF1 = 12 * log2(tfreqs(1, m) / tfreqs(1, 1));
                        p.DF2 = 12 * log2(tfreqs(2, m) / tfreqs(2, 1));
                        if (keepF > 2)
                            p.DF3 = 12 * log2(tfreqs(3, m) / tfreqs(3, 1));
                        else
                            p.DF3 = 0;
                        end
                        
                        %fprintf(logid, '%s ; %s ; %s ; %2.0f ; %2.0f ; %2.0f ; %2.0f ; %2.0f ; %2.0f ; %2.0f, %2.0f ; %2.0f\n', output, stimseq{1}, stimseq{2}, m-1, tfreqs(1, 1), tfreqs(1, end), tfreqs(1, m), p.DF1, tfreqs(2, 1), tfreqs(2, end), tfreqs(2, m), p.DF2);
                        if (keepF == 2)
                            savelog = sprintf('%s ; %s ; %s ; %2.0f ; %2.0f ; %2.0f ; %2.0f ; %2.3f ; %2.3f ; %2.0f ; %2.0f, %2.0f ; %2.3f ; %2.3f\n', output, stimstart, stimend, m-1, tfreqs(1, 1), tfreqs(1, end), tfreqs(1, m), coords(1, m), p.DF1, tfreqs(2, 1), tfreqs(2, end), tfreqs(2, m), coords(2, m), p.DF2);
                        else
                            savelog = sprintf('%s ; %s ; %s ; %2.0f ; %2.0f ; %2.0f ; %2.0f ; %2.3f ; %2.3f ; %2.0f ; %2.0f, %2.0f ; %2.3f ; %2.3f\n', output, stimstart, stimend, m-1, tfreqs(1, 1), tfreqs(1, end), tfreqs(1, m), coords(1, m), p.DF1, tfreqs(2, 1), tfreqs(2, end), tfreqs(2, m), coords(2, m), p.DF2, tfreqs(3, 1), tfreqs(3, end), tfreqs(3, m), coords(3, m), p.DF3);
                        end

                        fprintf(logid, '%s', savelog);
                        % + save graphs for future use
                        
                        if debug
                            % Save important information to file with
                            % filename, endpoints, freqs, distances...
                            tfreqs;
                            p.DF1;
                            p.DF2;
                            p.DF3;
                        end
                        
                        % manipulate_spectral_envelope.m applies on either a sound file or a
                        % corresponding vector. As it calls praat for formant extraction, it has to
                        % create a file from the vector.
                        
                        % Verify necessity to estimate formant freq + use extended options (+pb maxNformants)?
                        % should we locate the vowel mid-position for formant measurement?
                        % if so, need for list of pairs + positions;
                        % anyway need list of pairs for generating continua;
                        
                        % There's a weird issue with time output. Various options generate a
                        % different time vector
                        %[fx, r] = get_formants(audiopath, 'burg', 0.005, 5, 5500, 0.1);
                        %fx = get_formants(audiopath);
                        
                        if debug
                            [y, fs] = manipulate_spectral_envelope(sound, fs, p, true, true);
                        else
                            [y, fs] = manipulate_spectral_envelope(sound, fs, p, true, false);
                        end
                        
                        if wavestorage
                            outfile = [savedir, sprintf('%02i', i), '_', output]
                            audiowrite(outfile, y, fs);
                        end
                        %[y, fs] = manipulate_spectral_envelope(sound, fs, p);
                        
                        if debug
                            %s
                            %i
                            %j
                            %k
                            %l
                            %m
                            figure(226);
                            spectime(outfile, fparams);
                            pause(200/100)
                            snd = audioplayer(y, fs); 
                            play(snd);
                            pause(80/100);
                        end
                    end
                %end
             %end
            end
        end
   end
end

rmpath(pyPraat_path);

tEl = toc;
tEl/60

