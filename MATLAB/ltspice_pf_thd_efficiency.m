%% LTspice PF, THD and efficiency calculation
% Works with a tab/whitespace-delimited LTspice text export such as:
% time  V(ac_l,ac_n)  V(out,COM)  I(R_LOAD)  I(VAC)
%
% No Signal Processing Toolbox is required.

clear; clc;

%% ---------------- USER SETTINGS ----------------
fileName = 'Flyback_12v_9W_AC_corrected.txt';

% Exact LTspice trace names in the first row of the exported text file.
col.time = 'time';
col.vin  = 'V(ac_l,ac_n)';
col.iin  = 'I(VAC)';
col.vout = 'V(out,COM)';
col.iout = 'I(R_LOAD)';

fLine = 50;          % AC line frequency [Hz]
tStart = 0.460;      % analysis-window start [s]
tStop  = 0.500;      % analysis-window stop [s]
N_harmonics = 40;    % THD from harmonics 2...N (40 is common for line-current work)
makePlots = true;

% LTspice sign convention:
% I(VAC) is positive INTO the positive terminal of the voltage source.
% A source delivering power therefore usually has negative I(VAC).
inputCurrentSign  = -1;  % i_in drawn by converter = -I(VAC)
outputCurrentSign = +1;  % R_LOAD is from out to COM, so I(R_LOAD) is load current

%% ---------------- IMPORT LTSPICE TEXT FILE ----------------
[names, raw] = readLtspiceText(fileName);

idxTime = findColumn(names, col.time);
idxVin  = findColumn(names, col.vin);
idxIin  = findColumn(names, col.iin);
idxVout = findColumn(names, col.vout);
idxIout = findColumn(names, col.iout);

t    = raw(:, idxTime);
vin  = raw(:, idxVin);
iin  = inputCurrentSign  * raw(:, idxIin);
vout = raw(:, idxVout);
iout = outputCurrentSign * raw(:, idxIout);

% Remove invalid rows, sort by time, and remove duplicate time samples.
valid = all(isfinite([t, vin, iin, vout, iout]), 2);
t = t(valid);
Y = [vin(valid), iin(valid), vout(valid), iout(valid)];
[t, order] = sort(t);
Y = Y(order, :);
[t, uniqueIdx] = unique(t, 'stable');
Y = Y(uniqueIdx, :);

if tStart < t(1) || tStop > t(end) || tStop <= tStart
    error('Requested analysis window [%g, %g] s is outside the exported time range [%g, %g] s.', ...
          tStart, tStop, t(1), t(end));
end

% Force the window to contain exact tStart and tStop endpoints by interpolation.
[tw, Yw] = exactTimeWindow(t, Y, tStart, tStop);
vin  = Yw(:,1);
iin  = Yw(:,2);
vout = Yw(:,3);
iout = Yw(:,4);

nCycles = (tStop - tStart) * fLine;
if abs(nCycles - round(nCycles)) > 1e-8
    warning(['The selected window is %.9g line cycles, not an integer number. ' ...
             'Choose tStop-tStart = integer/fLine to avoid spectral leakage.'], nCycles);
end

%% ---------------- TIME-DOMAIN POWER AND RMS ----------------
Vin_rms = timeRms(tw, vin);
Iin_rms = timeRms(tw, iin);
Vout_rms = timeRms(tw, vout);
Iout_rms = timeRms(tw, iout);

Vin_avg = timeAverage(tw, vin);
Iin_avg = timeAverage(tw, iin);
Vout_avg = timeAverage(tw, vout);
Iout_avg = timeAverage(tw, iout);

Pin  = timeAverage(tw, vin .* iin);   % true AC input power [W]
Pout = timeAverage(tw, vout .* iout); % exact average load power [W]
Sin  = Vin_rms * Iin_rms;             % apparent input power [VA]
truePF = Pin / Sin;
efficiency_pct = 100 * Pout / Pin;

%% ---------------- HARMONIC ANALYSIS ----------------
% Harmonic coefficients are calculated by direct time integration, so the
% original nonuniform LTspice time samples do not need to be resampled.
Hvin = harmonicAnalysis(tw, vin, fLine, N_harmonics);
Hiin = harmonicAnalysis(tw, iin, fLine, N_harmonics);

% Harmonic THD: only integer line harmonics 2...N_harmonics.
V_THD_harm_pct = 100 * sqrt(sum(Hvin.rms(2:end).^2)) / Hvin.rms(1);
I_THD_harm_pct = 100 * sqrt(sum(Hiin.rms(2:end).^2)) / Hiin.rms(1);

% Residual THD: all non-DC, non-fundamental content present in the exported
% waveform. This also includes switching ripple and numerical noise.
V_THD_residual_pct = 100 * residualTHD(Vin_rms, Hvin.dc, Hvin.rms(1));
I_THD_residual_pct = 100 * residualTHD(Iin_rms, Hiin.dc, Hiin.rms(1));

% Fundamental displacement power factor. Unlike true PF, this only measures
% the phase angle between the 50-Hz voltage and current fundamentals.
V1 = Hvin.phasorRms(1);
I1 = Hiin.phasorRms(1);
P1 = real(V1 * conj(I1));
displacementPF = P1 / (abs(V1) * abs(I1));
phaseDifference_deg = angle(V1 * conj(I1)) * 180/pi;
currentDistortionFactor = abs(I1) / Iin_rms;

%% ---------------- STEADY-STATE CHECK ----------------
voutDrift = vout(end) - vout(1);
voutDrift_pct = 100 * voutDrift / max(abs(Vout_avg), eps);

%% ---------------- DISPLAY RESULTS ----------------
fprintf('\nAnalysis window: %.6f s to %.6f s (%.3f line cycles)\n', ...
        tStart, tStop, nCycles);
fprintf('------------------------------------------------------------\n');
fprintf('Input voltage RMS                 = %12.6g V\n', Vin_rms);
fprintf('Input current RMS                 = %12.6g A\n', Iin_rms);
fprintf('Average input power               = %12.6g W\n', Pin);
fprintf('Input apparent power              = %12.6g VA\n', Sin);
fprintf('TRUE power factor                 = %12.6f\n', truePF);
fprintf('Fundamental displacement PF       = %12.6f\n', displacementPF);
fprintf('V-I fundamental phase difference  = %12.6f deg\n', phaseDifference_deg);
fprintf('Current distortion factor I1/Irms = %12.6f\n', currentDistortionFactor);
fprintf('\n');
fprintf('Input-current THD, harmonics 2-%d = %12.6f %%\n', ...
        N_harmonics, I_THD_harm_pct);
fprintf('Input-current residual THD         = %12.6f %%\n', ...
        I_THD_residual_pct);
fprintf('Input-voltage THD, harmonics 2-%d = %12.6f %%\n', ...
        N_harmonics, V_THD_harm_pct);
fprintf('Input-voltage residual THD         = %12.6f %%\n', ...
        V_THD_residual_pct);
fprintf('\n');
fprintf('Average output voltage             = %12.6g V\n', Vout_avg);
fprintf('Average output current             = %12.6g A\n', Iout_avg);
fprintf('Average output power               = %12.6g W\n', Pout);
fprintf('EFFICIENCY                          = %12.6f %%\n', efficiency_pct);
fprintf('Output-voltage drift over window   = %12.6g V (%g %%)\n', ...
        voutDrift, voutDrift_pct);

if Pin <= 0
    warning(['Average input power is not positive. Check inputCurrentSign ' ...
             'and the voltage/current reference directions.']);
end
if Pout <= 0
    warning(['Average output power is not positive. Check outputCurrentSign ' ...
             'and the load current reference direction.']);
end
if abs(voutDrift_pct) > 0.5
    warning(['Output voltage changes by more than 0.5%% during the analysis window. ' ...
             'Extend the simulation and measure after steady state is reached.']);
end

%% Results are also left in convenient MATLAB variables.
results = struct(...
    'Vin_rms_V', Vin_rms, ...
    'Iin_rms_A', Iin_rms, ...
    'Pin_W', Pin, ...
    'Sin_VA', Sin, ...
    'TruePF', truePF, ...
    'DisplacementPF', displacementPF, ...
    'PhaseDifference_deg', phaseDifference_deg, ...
    'InputCurrentTHD_harmonic_pct', I_THD_harm_pct, ...
    'InputCurrentTHD_residual_pct', I_THD_residual_pct, ...
    'InputVoltageTHD_harmonic_pct', V_THD_harm_pct, ...
    'Pout_W', Pout, ...
    'Efficiency_pct', efficiency_pct, ...
    'Vout_avg_V', Vout_avg, ...
    'Iout_avg_A', Iout_avg, ...
    'VoutDrift_pct', voutDrift_pct);

harmonicNumber = (1:N_harmonics).';
harmonics = table(harmonicNumber, Hvin.rms, 100*Hvin.rms/Hvin.rms(1), ...
                  Hiin.rms, 100*Hiin.rms/Hiin.rms(1), ...
    'VariableNames', {'Harmonic','V_RMS_V','V_pct_of_fundamental', ...
                      'I_RMS_A','I_pct_of_fundamental'});

disp(' ');
disp('First 20 harmonic rows:');
disp(harmonics(1:min(20,height(harmonics)), :));

% Uncomment to save all harmonic values:
% writetable(harmonics, 'ltspice_harmonics.csv');

%% ---------------- OPTIONAL PLOTS ----------------
if makePlots
    figure('Name','LTspice input waveforms');
    yyaxis left
    plot(tw, vin, 'LineWidth', 1.0);
    ylabel('Input voltage [V]');
    yyaxis right
    plot(tw, iin, 'LineWidth', 1.0);
    ylabel('Input current [A]');
    xlabel('Time [s]');
    title(sprintf('Input voltage and current, PF = %.4f', truePF));
    grid on;

    figure('Name','Input-current harmonics');
    stem(harmonicNumber, 100*Hiin.rms/Hiin.rms(1), 'filled');
    xlabel('Harmonic number');
    ylabel('Current harmonic [% of fundamental]');
    title(sprintf('Input-current harmonics; THD_{2-%d} = %.2f%%', ...
                  N_harmonics, I_THD_harm_pct));
    grid on;
end

%% ================= LOCAL FUNCTIONS =================
function [names, data] = readLtspiceText(fileName)
    fid = fopen(fileName, 'r');
    if fid < 0
        error('Cannot open file: %s', fileName);
    end
    cleaner = onCleanup(@() fclose(fid));

    header = fgetl(fid);
    if ~ischar(header)
        error('File is empty: %s', fileName);
    end

    % LTspice normally exports tab-separated text. Whitespace fallback is
    % included for files that were edited or copied through another program.
    if contains(header, sprintf('\t'))
        names = strtrim(strsplit(strtrim(header), sprintf('\t')));
    else
        names = regexp(strtrim(header), '\s+', 'split');
    end

    nColumns = numel(names);
    formatSpec = repmat('%f', 1, nColumns);
    C = textscan(fid, formatSpec, ...
        'Delimiter', sprintf('\t '), ...
        'MultipleDelimsAsOne', true, ...
        'CollectOutput', true, ...
        'ReturnOnError', false);

    data = C{1};
    if isempty(data) || size(data,2) ~= nColumns
        error(['Could not parse the numeric data. Found %d header columns but %d ' ...
               'numeric columns. Make sure the file is a normal LTspice text export.'], ...
              nColumns, size(data,2));
    end
end

function idx = findColumn(names, requestedName)
    idx = find(strcmpi(strtrim(names), strtrim(requestedName)), 1);
    if isempty(idx)
        fprintf('Available columns in the file:\n');
        for k = 1:numel(names)
            fprintf('  %2d: %s\n', k, names{k});
        end
        error('Column "%s" was not found. Copy its exact name into the settings section.', ...
              requestedName);
    end
end

function [tw, Yw] = exactTimeWindow(t, Y, tStart, tStop)
    inside = t > tStart & t < tStop;
    tw = [tStart; t(inside); tStop];
    Ystart = interp1(t, Y, tStart, 'linear');
    Ystop  = interp1(t, Y, tStop,  'linear');
    Yw = [Ystart; Y(inside,:); Ystop];
end

function value = timeAverage(t, x)
    value = trapz(t, x) / (t(end) - t(1));
end

function value = timeRms(t, x)
    value = sqrt(max(trapz(t, x.^2) / (t(end) - t(1)), 0));
end

function H = harmonicAnalysis(t, x, f0, nMax)
    duration = t(end) - t(1);
    H.dc = trapz(t, x) / duration;
    H.rms = zeros(nMax,1);
    H.phase_rad = zeros(nMax,1);
    H.phasorRms = complex(zeros(nMax,1));

    for n = 1:nMax
        omega = 2*pi*n*f0;
        aCos = (2/duration) * trapz(t, x .* cos(omega*t));
        bSin = (2/duration) * trapz(t, x .* sin(omega*t));

        % x_n(t) = aCos*cos(wt) + bSin*sin(wt)
        % RMS phasor convention: X = (aCos - j*bSin)/sqrt(2)
        H.phasorRms(n) = (aCos - 1j*bSin) / sqrt(2);
        H.rms(n) = abs(H.phasorRms(n));
        H.phase_rad(n) = angle(H.phasorRms(n));
    end
end

function thdRatio = residualTHD(totalRms, dcValue, fundamentalRms)
    nonFundamentalSquared = totalRms^2 - dcValue^2 - fundamentalRms^2;
    thdRatio = sqrt(max(nonFundamentalSquared, 0)) / max(fundamentalRms, eps);
end
