function burst_properties = threshold_high_amplitude_events(time, signal, freq, params)

% THRESHOLD_HIGH_AMPLITUDE_EVENTS Detects and characterizes high-amplitude oscillatory bursts.
%
%   burst_properties = THRESHOLD_HIGH_AMPLITUDE_EVENTS(time, signal, freq, params)
%
%   This function detects transient high-amplitude bursts in a 1D signal
%   within a specified frequency band and returns an output structure with
%   burst timing, durations, and summary measures. The signal is bandpass-
%   filtered, the analytic amplitude is obtained via the Hilbert transform,
%   and bursts are defined as contiguous periods where the amplitude exceeds
%   a user-defined percentile threshold. Bursts are then filtered for
%   minimum duration and distance from the signal edges, with optional
%   merging of temporally adjacent events and visualization of burst
%   statistics. Conceptually similar approaches using Hilbert amplitude and
%   thresholding are common in neural oscillation burst detection. [web:18][web:19]
%
%   INPUTS:
%       time    - 1D vector of time points in seconds, uniformly sampled,
%                 same length as signal.
%       signal  - 1D vector of the time-domain signal.
%       freq    - Two-element vector [f_low, f_high] specifying the
%                 bandpass range in Hz used for burst detection.
%       params  - (Optional) struct with fields:
%           .thresh       : Percentile for amplitude threshold (default: 75),
%                           applied to the Hilbert envelope.
%           .edge_space   : Edge exclusion window in seconds; bursts whose
%                           samples fall within this distance of the start
%                           or end of the signal are discarded (default: 1).
%           .edge_cut     : Duration in seconds to trim from both ends of
%                           the signal before bandpass filtering, to reduce
%                           edge artifacts (default: 0).
%           .merge        : Logical flag indicating whether to merge bursts
%                           whose temporal separation is shorter than the
%                           minimum burst duration (default: false).
%           .burst_epoch  : String specifying which burst time point to use
%                           as the burst index:
%                              'first_crossing'  - first above-threshold sample
%                              'max_amplitude'   - sample with maximal absolute amplitude
%                              'end_crossing'    - last above-threshold sample
%                           (default: 'first_crossing').
%           .vis          : Logical flag; if true, plots histograms of burst
%                           durations (pre/post filtering) and the final
%                           burst vector over time (default: true).
%
%   OUTPUT:
%       burst_properties  - Struct containing:
%           .raw_burst_vector      : Logical vector (length of filtered
%                                   segment) indicating above-threshold
%                                   samples before duration/edge filtering.
%           .detection_threshold   : Scalar amplitude threshold (same units
%                                   as Hilbert envelope) corresponding to
%                                   the specified percentile.
%           .raw_CCObj             : bwconncomp output for initial burst
%                                   segments before merging or filtering.
%           .merged_CCObj          : (If merge == true) bwconncomp output
%                                   after merging adjacent events.
%           .postproc_burst_vector : Logical vector (same length as filtered
%                                   segment) marking samples belonging to
%                                   bursts that pass duration and edge
%                                   criteria.
%           .postproc_event_idx    : Vector of global sample indices in the
%                                   original signal corresponding to each
%                                   burst’s representative time point
%                                   (defined by burst_epoch).
%           .burst_proportion      : Fraction of samples classified as
%                                   bursts (fractional occupancy).
%           .burst_rate            : Burst rate in bursts per second.
%           .burst_durations       : Vector of accepted burst durations in
%                                   milliseconds.
%           .used_params           : Struct echoing all parameters actually
%                                   used (including derived fields such as
%                                   sampling frequency and minimum burst
%                                   duration).
%
%   ALGORITHM (OVERVIEW):
%       1) Estimate sampling frequency from the time vector.
%       2) Optionally trim edges and bandpass-filter the signal in [f_low, f_high].
%       3) Compute the Hilbert envelope and derive a percentile-based
%          amplitude threshold.
%       4) Threshold the envelope to obtain a binary burst vector.
%       5) Optionally merge temporally adjacent bursts whose separation is
%          shorter than the minimum burst duration.
%       6) Identify contiguous burst segments and discard those that:
%            - Are shorter than one cycle of the lowest frequency in freq.
%            - Occur too close to the signal edges (edge_space).
%       7) For accepted bursts, compute durations and a representative index
%          (first sample, peak, or last sample, depending on burst_epoch).
%       8) Compute burst proportion and burst rate and optionally visualize
%          duration distributions and burst time course.
%
%   REQUIREMENTS:
%       - Assumes uniformly sampled data.
%       - Uses functions such as BANDPASS, HILBERT, BWCONNCOMP, HISTOGRAM,
%         and PLOT from MATLAB and the Image Processing Toolbox where
%         applicable.
%
%   See also: bandpass, hilbert, bwconncomp

%% default params
thresh = 75; %
edge_space = 1;
edge_cut = 0;
merge = false;
burst_epoch = 'first_crossing';
vis = true;
%% parse params
if exist('params', 'var') && ~isempty(params)

    if isfield(params, 'thresh')
        thresh = params.thresh;
    end

    if isfield(params, 'vis')
        vis = params.vis;
    end

    if isfield(params, 'edge_space')
        edge_space = params.edge_space;
    end

    if isfield(params, 'edge_cut')
        edge_cut = params.edge_cut;
    end

    if isfield(params, 'burst_epoch')
        burst_epoch = params.burst_epoch;
    end

    if isfield(params, 'merge')
        merge = params.merge;
    end
end
% Set up output structure
burst_properties = struct();

%% extract high amplitude events
% get sampling frequency
sfreq = round(1/(time(2)-time(1)));
edge_space_samples = edge_space*sfreq;

%% bandpass filter the signal
filtered_signal = bandpass(signal(sfreq*edge_cut+1:end-sfreq*edge_cut), freq, sfreq);

%% threshold signal
% Amplitude as a default for AG Florin
% Here one can also take the square of the absolute to get the instantaneous
% power
hilbert_env = abs(hilbert(filtered_signal));
% Get detection threshold. 
detection_thresh = prctile(hilbert_env, thresh);

% Apply detection threshold to amplitude timeseires -> Amplitude bigger
% then threshold -> burst.
burst_vec = hilbert_env >= detection_thresh;
burst_properties.raw_burst_vector = burst_vec;
burst_properties.detection_threshold = detection_thresh;
%% check for validity
%check for minimum duration of at least one cycle of the lowest frequency
min_burst_duration = (1/min(freq)/(1/sfreq));
CC = bwconncomp(burst_vec);
burst_properties.raw_CCObj = CC;

% Merge events if specified
if merge
    for ii=1:length(CC.PixelIdxList)
        if ii == 1
            continue
        end
        prev_end_idx = max(CC.PixelIdxList{ii-1});
        curr_start_idx = min(CC.PixelIdxList{ii});
        % If start of current event minus the minimum burst duration
        % overlaps with the end of the prior event combine them to one
        % burst
        if curr_start_idx - min_burst_duration < prev_end_idx
            burst_vec(min(CC.PixelIdxList{ii-1}):max(CC.PixelIdxList{ii})) = true;
        end
    end
    CC = bwconncomp(burst_vec);
    burst_properties.merged_CCObj = CC;
end

event_length_pre = zeros(length(CC.NumObjects), 1);
event_length_post = [];
event_idx = [];

for i = 1:CC.NumObjects
    idx = CC.PixelIdxList{i};
    %get burst length in ms before filtering too short ones out
    event_length_pre(i) = (numel(idx)/sfreq)*1000;
    
    % make sure burst are of sufficient length and don't appear to close to
    % the edge
    if numel(idx) < min_burst_duration || min(idx) < edge_space_samples || max(idx) > length(signal)-edge_space_samples 
        burst_vec(idx) = 0;
    else
        %get length of suffiently long bursts
        event_length_post(end+1) = (numel(idx)/sfreq)*1000; %converting duration of burst to ms
        
        event = signal(idx+sfreq*edge_cut);
        event = event - mean(event);
        
        % Get burst idx depending on specified method
        if strcmp(burst_epoch, 'first_crossing')
            event_idx(end+1) = idx(1); % Take beginning of the burst
        elseif strcmp(burst_epoch, 'max_amplitude')
            [~, relative_idx] = max(abs(event)); % find peak burst activity
            event_idx(end+1) = idx(relative_idx); % get global index of that peak idx.
        elseif strcmp(burst_epoch, 'end_crossing')
            event_idx(end+1) = idx(end); % Take end of the burst
        end

        %remove nan's from array
        event_idx = event_idx(~isnan(event_idx));
     end
end

% Calculate burst descriptives
burst_proportion = sum(burst_vec) / length(burst_vec); % Same as "fractional occupancy"
burst_rate = length(event_idx)/length(signal)/sfreq; % bursts per second

% Fill output structure
burst_properties.postproc_burst_vector = burst_vec;
burst_properties.postproc_event_idx = event_idx;
burst_properties.burst_proportion = burst_proportion;
burst_properties.burst_rate = burst_rate;
burst_properties.burst_durations = event_length_post; % already in ms

used_params = struct( ...
    'thresh', thresh, ...
    'edge_space', edge_space, ...
    'edge_cut', edge_cut, ...
    'merge', merge, ...
    'burst_epoch', burst_epoch, ...
    'vis', vis, ...
    'sfreq', sfreq, ...
    'min_burst_duration', min_burst_duration);

burst_properties.used_params = used_params;


%% visualize distribution of events length
if vis && length(event_idx) > 1
    f = figure();
    tiledlayout(3, 1);
    fig_width = 1200;
    fig_height = 1400;
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');
    set(gcf, 'Position',  [1, 1, 100 + fig_width, 100 + fig_height]);   % Resize fig window

    % First tile: unfiltered burst length histogram
    nexttile();
    hold("on");
    histogram(event_length_pre, floor(3.49*std(event_length_pre)*length(event_length_pre)^(1/3)), ...
        'BinWidth', 3);
    title(sprintf('Burst length distribution without filtering out too short ones (n=%i).', length(event_length_pre)));

    nexttile();
    hold("on");
    histogram(event_length_post, floor(3.49*std(event_length_post)*length(event_length_post)^(1/3)), ...
        'BinWidth', 3);
    title(sprintf('Burst length distribution with filtering out too short and overlapping ones (n=%i).', length(event_idx)));
    
    nexttile();
    hold("on")
    plot(time, burst_vec);
    title("Burst distribution (filtered) over time.");
end
