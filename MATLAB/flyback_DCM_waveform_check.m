%% Standalone DCM verification plot for LTspice flyback
% Export two LTspice text files over the SAME short time interval:
%
% dcm_currents.txt:
%   time   I(Lp)   I(Ls)   I(D§Schottky)
%
% dcm_switching.txt:
%   time   V(Drain,CS)   V(PWM)
%
% This script plots the waveforms needed to visually demonstrate DCM.
% DCM evidence:
%   1) PWM turns OFF.
%   2) Secondary/diode current conducts and decays.
%   3) Secondary/diode current reaches ~0 A.
%   4) It remains ~0 A for a finite interval.
%   5) Only then does the next PWM turn-on occur.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
currentFile   = 'Flyback_12v_9W_AC_corrected_dcm_currents.txt';
switchingFile = 'Flyback_12v_9W_AC_corrected_dcm_switching.txt';

% Exact LTspice trace names
col.time = 'time';
col.ip   = 'I(Lp)';
col.is   = 'I(Ls)';
col.id   = 'I(Schottky)';
col.vds  = 'V(Drain,CS)';
col.pwm  = 'V(PWM)';

% Thresholds used only to mark the end of demagnetization.
% Adjust if LTspice ringing/noise prevents the currents from reaching
% numerically exact zero.
currentThreshold_A = 1e-3;    % 1 mA
pwmThreshold_V     = 0.5;     % PWM considered ON above this level

%% ---------------- IMPORT ----------------
[names1, data1] = readLtspiceText(currentFile);
[names2, data2] = readLtspiceText(switchingFile);

t1 = data1(:, findColumn(names1, col.time));
ip = data1(:, findColumn(names1, col.ip));
isec = data1(:, findColumn(names1, col.is));
id = data1(:, findColumn(names1, col.id));

t2 = data2(:, findColumn(names2, col.time));
vds = data2(:, findColumn(names2, col.vds));
pwm = data2(:, findColumn(names2, col.pwm));

% Remove invalid rows
ok1 = all(isfinite([t1 ip isec id]), 2);
t1 = t1(ok1);
ip = ip(ok1);
isec = isec(ok1);
id = id(ok1);

ok2 = all(isfinite([t2 vds pwm]), 2);
t2 = t2(ok2);
vds = vds(ok2);
pwm = pwm(ok2);

% Sort and remove duplicate time samples
[t1, ia] = unique(t1, 'sorted');
ip = ip(ia);
isec = isec(ia);
id = id(ia);

[t2, ib] = unique(t2, 'sorted');
vds = vds(ib);
pwm = pwm(ib);

% Interpolate switching signals onto the current-file time vector.
% Both files should have been exported over the same LTspice time window.
if t1(1) < t2(1) || t1(end) > t2(end)
    error(['The two DCM files do not cover the same time range. ' ...
           'Export them over the same interval.']);
end

vds_i = interp1(t2, vds, t1, 'linear');
pwm_i = interp1(t2, pwm, t1, 'previous');

t = t1;

%% ---------------- FIND PWM TURN-ON EVENTS ----------------
pwmOn = pwm_i > pwmThreshold_V;
turnOnIdx = find(diff(pwmOn) == 1) + 1;

if numel(turnOnIdx) < 2
    error('Could not find at least two PWM turn-on events.');
end

%% ---------------- FIND A REPRESENTATIVE DCM CYCLE ----------------
% For each switching period, find the last point where either secondary
% current or diode current is above the threshold. The interval from there
% to the next turn-on is the DCM zero-current interval.

bestCycle = NaN;
bestDeadTime = -Inf;
demagIdxBest = NaN;

for k = 1:numel(turnOnIdx)-1
    i0 = turnOnIdx(k);
    i1 = turnOnIdx(k+1);

    active = abs(isec(i0:i1)) > currentThreshold_A | ...
             abs(id(i0:i1))   > currentThreshold_A;

    lastActiveRel = find(active, 1, 'last');

    if isempty(lastActiveRel)
        continue;
    end

    demagIdx = i0 + lastActiveRel - 1;

    if demagIdx < i1
        deadTime = t(i1) - t(demagIdx);

        if deadTime > bestDeadTime
            bestDeadTime = deadTime;
            bestCycle = k;
            demagIdxBest = demagIdx;
        end
    end
end

if isnan(bestCycle)
    error(['No zero-current interval was detected. Increase the exported ' ...
           'window or adjust currentThreshold_A.']);
end

nextTurnOnIdx = turnOnIdx(bestCycle + 1);

%% ---------------- DISPLAY RESULT ----------------
fprintf('\nDCM waveform check\n');
fprintf('--------------------------------------------------\n');
fprintf('Detected PWM periods                 = %d\n', numel(turnOnIdx)-1);
fprintf('Displayed zero-current interval      = %.3f us\n', ...
        1e6*(t(nextTurnOnIdx)-t(demagIdxBest)));
fprintf('Secondary current at next turn-on    = %.6g A\n', isec(nextTurnOnIdx));
fprintf('Diode current at next turn-on        = %.6g A\n', id(nextTurnOnIdx));
fprintf('\n');
fprintf(['DCM evidence: secondary/diode current reaches approximately zero\n' ...
         'before the next PWM turn-on and remains there for a finite interval.\n']);

%% ---------------- PLOT A FEW SWITCHING CYCLES ----------------
% Show one cycle before and one cycle after the chosen interval if possible.
leftCycle  = max(bestCycle-1, 1);
rightCycle = min(bestCycle+2, numel(turnOnIdx));

plotStart = turnOnIdx(leftCycle);
plotStop  = turnOnIdx(rightCycle);

tp = 1e6*(t(plotStart:plotStop) - t(plotStart));  % microseconds
idx = plotStart:plotStop;

demag_us = 1e6*(t(demagIdxBest) - t(plotStart));
nextOn_us = 1e6*(t(nextTurnOnIdx) - t(plotStart));

figure('Name','Flyback DCM verification','Color','w');

subplot(5,1,1)
plot(tp, ip(idx), 'LineWidth', 1.0);
ylabel('I_L_p [A]');
grid on;
title('Flyback DCM verification');

subplot(5,1,2)
plot(tp, isec(idx), 'LineWidth', 1.0);
ylabel('I_L_s [A]');
grid on;

subplot(5,1,3)
plot(tp, id(idx), 'LineWidth', 1.0);
ylabel('I_D [A]');
grid on;

subplot(5,1,4)
plot(tp, vds_i(idx), 'LineWidth', 1.0);
ylabel('V_{DS} [V]');
grid on;

subplot(5,1,5)
plot(tp, pwm_i(idx), 'LineWidth', 1.0);
ylabel('PWM [V]');
xlabel('Time [\mus]');
grid on;

% Mark end of energy transfer and next PWM turn-on on every subplot.
for s = 1:5
    subplot(5,1,s)
    xline(demag_us, '--', 'demag end', ...
        'LabelVerticalAlignment','bottom');
    xline(nextOn_us, ':', 'next turn-on', ...
        'LabelVerticalAlignment','bottom');
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
        error(['Could not parse %s. Found %d header columns but %d numeric columns.'], ...
              fileName, nColumns, size(data,2));
    end
end

function idx = findColumn(names, requestedName)
    idx = find(strcmpi(strtrim(names), strtrim(requestedName)), 1);

    if isempty(idx)
        fprintf('Available columns:\n');
        for k = 1:numel(names)
            fprintf('  %2d: %s\n', k, names{k});
        end
        error('Column "%s" was not found.', requestedName);
    end
end
