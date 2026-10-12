%% Flyback converter: efficiency, PF and THD with a switched load
% Required LTspice export columns:
%   time  V(ac_l,ac_n)  I(VAC)  V(out,COM)  I(R_LOAD)  I(R_VARYING)
%
% Optional column:
%   V(N004,COM)
%
% If V(N004,COM) is present, the varying-resistor power is calculated as
% V(N004,COM)*I(R_VARYING). Otherwise it is calculated as
% I(R_VARYING)^2*Rvarying.
%
% No Signal Processing Toolbox is required.

clear; clc;

%% ---------------- USER SETTINGS ----------------
fileName = 'Flyback_12v_9W_AC_PFC_TIME-VARYING_17.3_5.76.txt';

% Exact trace names in the first row of the LTspice text export.
col.time    = 'time';
col.vin     = 'V(ac_l,ac_n)';
col.iin     = 'I(VAC)';
col.vout    = 'V(out,COM)';
col.iFixed  = 'I(R_LOAD)';
col.iVary   = 'I(R_VARYING)';
col.vVary   = 'V(N004,COM)';       % optional

fLine = 50;                        % line frequency [Hz]
tStart = 0.460;                    % analysis start [s]
tStop  = 0.500;                    % analysis end [s]
N_harmonics = 40;                  % current harmonics 2...40
Rvarying = 5.76;                   % used only if col.vVary is absent [ohm]
makePlots = true;

% LTspice defines I(VAC) as current entering the source's positive terminal.
% Therefore the current drawn by the converter is normally -I(VAC).
inputCurrentSign = -1;

%% ---------------- IMPORT LTSPICE DATA ----------------
[names, raw] = readLtspiceText(fileName);

idx.time   = requireColumn(names, col.time);
idx.vin    = requireColumn(names, col.vin);
idx.iin    = requireColumn(names, col.iin);
idx.vout   = requireColumn(names, col.vout);
idx.iFixed = requireColumn(names, col.iFixed);
idx.iVary  = requireColumn(names, col.iVary);
idx.vVary  = optionalColumn(names, col.vVary);

t      = raw(:,idx.time);
vin    = raw(:,idx.vin);
iin    = inputCurrentSign*raw(:,idx.iin);
vout   = raw(:,idx.vout);
iFixed = raw(:,idx.iFixed);
iVary  = raw(:,idx.iVary);

hasVaryingVoltage = ~isempty(idx.vVary);
if hasVaryingVoltage
    vVary = raw(:,idx.vVary);
    allSignals = [vin, iin, vout, iFixed, iVary, vVary];
else
    allSignals = [vin, iin, vout, iFixed, iVary];
end

% Remove invalid rows, sort by time and remove duplicate timestamps.
valid = isfinite(t) & all(isfinite(allSignals),2);
t = t(valid);
allSignals = allSignals(valid,:);
[t,order] = sort(t);
allSignals = allSignals(order,:);
[t,uniqueRows] = unique(t,'stable');
allSignals = allSignals(uniqueRows,:);

if tStop <= tStart || tStart < t(1) || tStop > t(end)
    error(['Analysis window [%g, %g] s is outside the exported range ' ...
           '[%g, %g] s.'],tStart,tStop,t(1),t(end));
end

% Interpolate exact start and stop points so every integral uses precisely
% the requested window, even when LTspice used nonuniform timesteps.
[tw,Yw] = exactTimeWindow(t,allSignals,tStart,tStop);
vin    = Yw(:,1);
iin    = Yw(:,2);
vout   = Yw(:,3);
iFixed = Yw(:,4);
iVary  = Yw(:,5);
if hasVaryingVoltage
    vVary = Yw(:,6);
end

duration = tStop-tStart;
nLineCycles = duration*fLine;
if abs(nLineCycles-round(nLineCycles)) > 1e-8
    warning(['The analysis interval contains %.9g line cycles. Use an ' ...
             'integer number of cycles for reliable PF and THD.'],nLineCycles);
end

%% ---------------- INSTANTANEOUS AND AVERAGE POWER ----------------
pIn = vin.*iin;
pFixed = vout.*iFixed;

if hasVaryingVoltage
    pVary = vVary.*iVary;
    varyingPowerMethod = 'V(N004,COM) * I(R_VARYING)';
else
    % This is the power dissipated only in R_VARYING. It does not include
    % the small S_LOAD conduction loss.
    pVary = iVary.^2*Rvarying;
    varyingPowerMethod = sprintf('I(R_VARYING)^2 * %.6g ohm',Rvarying);
end

pOut = pFixed+pVary;

Ein       = trapz(tw,pIn);
Efixed    = trapz(tw,pFixed);
Evary     = trapz(tw,pVary);
Eout      = Efixed+Evary;
Pin       = Ein/duration;
Pfixed    = Efixed/duration;
Pvary     = Evary/duration;
Pout      = Eout/duration;
efficiency_pct = 100*Eout/Ein;

%% ---------------- RMS, PF AND HARMONICS ----------------
Vin_rms = timeRms(tw,vin);
Iin_rms = timeRms(tw,iin);
Vout_avg = timeAverage(tw,vout);
Vout_min = min(vout);
Vout_max = max(vout);
Sin = Vin_rms*Iin_rms;
truePF = Pin/Sin;

Hvin = harmonicAnalysis(tw,vin,fLine,N_harmonics);
Hiin = harmonicAnalysis(tw,iin,fLine,N_harmonics);

V_THD_harm_pct = 100*sqrt(sum(Hvin.rms(2:end).^2))/max(Hvin.rms(1),eps);
I_THD_harm_pct = 100*sqrt(sum(Hiin.rms(2:end).^2))/max(Hiin.rms(1),eps);
V_THD_residual_pct = 100*residualTHD(Vin_rms,Hvin.dc,Hvin.rms(1));
I_THD_residual_pct = 100*residualTHD(Iin_rms,Hiin.dc,Hiin.rms(1));

V1 = Hvin.phasorRms(1);
I1 = Hiin.phasorRms(1);
displacementPF = real(V1*conj(I1))/(max(abs(V1)*abs(I1),eps));
phaseDifference_deg = angle(V1*conj(I1))*180/pi;
currentDistortionFactor = abs(I1)/max(Iin_rms,eps);

%% ---------------- PERIODIC-STATE CHECK ----------------
voutChange = vout(end)-vout(1);
voutChange_pct = 100*voutChange/max(abs(Vout_avg),eps);

%% ---------------- DISPLAY RESULTS ----------------
fprintf('\nAnalysis window: %.6f to %.6f s (%.3f line cycles)\n', ...
        tStart,tStop,nLineCycles);
fprintf('------------------------------------------------------------\n');
fprintf('Input voltage RMS                   = %12.6g V\n',Vin_rms);
fprintf('Input current RMS                   = %12.6g A\n',Iin_rms);
fprintf('Average input power                 = %12.6g W\n',Pin);
fprintf('Input apparent power                = %12.6g VA\n',Sin);
fprintf('True power factor                   = %12.6f\n',truePF);
fprintf('Fundamental displacement PF         = %12.6f\n',displacementPF);
fprintf('V-I fundamental phase difference    = %12.6f deg\n',phaseDifference_deg);
fprintf('Current distortion factor I1/Irms   = %12.6f\n',currentDistortionFactor);
fprintf('\n');
fprintf('Input-current THD, harmonics 2-%d   = %12.6f %%\n', ...
        N_harmonics,I_THD_harm_pct);
fprintf('Input-current residual THD          = %12.6f %%\n',I_THD_residual_pct);
fprintf('Input-voltage THD, harmonics 2-%d   = %12.6f %%\n', ...
        N_harmonics,V_THD_harm_pct);
fprintf('Input-voltage residual THD          = %12.6f %%\n',V_THD_residual_pct);
fprintf('\n');
fprintf('Average output voltage              = %12.6g V\n',Vout_avg);
fprintf('Minimum / maximum output voltage    = %12.6g / %.6g V\n', ...
        Vout_min,Vout_max);
fprintf('Fixed-load power                    = %12.6g W\n',Pfixed);
fprintf('Varying-load power                  = %12.6g W\n',Pvary);
fprintf('Total useful output power           = %12.6g W\n',Pout);
fprintf('Efficiency                          = %12.6f %%\n',efficiency_pct);
fprintf('Varying-load calculation            = %s\n',varyingPowerMethod);
fprintf('Output-voltage endpoint change      = %12.6g V (%g %%)\n', ...
        voutChange,voutChange_pct);

if Pin <= 0
    warning(['Input power is not positive. Check inputCurrentSign and the ' ...
             'VAC voltage/current reference directions.']);
end
if Pfixed < 0 || Pvary < 0
    warning(['A load power is negative. Check the LTspice resistor current ' ...
             'directions and exported trace names.']);
end
if efficiency_pct > 100.5
    warning(['Efficiency exceeds 100%%. Check that the run matches the ' ...
             'netlist, the power signs, and whether stored or external-source ' ...
             'energy is contributing during the selected interval.']);
end
if abs(voutChange_pct) > 0.5
    warning(['Output voltage changes by more than 0.5%% between the interval ' ...
             'endpoints. Use a full periodic window or evaluate stored energy.']);
end

%% Results remain available in the MATLAB workspace.
results = struct( ...
    'tStart_s',tStart, ...
    'tStop_s',tStop, ...
    'Vin_rms_V',Vin_rms, ...
    'Iin_rms_A',Iin_rms, ...
    'Pin_W',Pin, ...
    'Sin_VA',Sin, ...
    'TruePF',truePF, ...
    'DisplacementPF',displacementPF, ...
    'PhaseDifference_deg',phaseDifference_deg, ...
    'InputCurrentTHD_harmonic_pct',I_THD_harm_pct, ...
    'InputCurrentTHD_residual_pct',I_THD_residual_pct, ...
    'InputVoltageTHD_harmonic_pct',V_THD_harm_pct, ...
    'Vout_avg_V',Vout_avg, ...
    'Pfixed_W',Pfixed, ...
    'Pvarying_W',Pvary, ...
    'Pout_total_W',Pout, ...
    'Efficiency_pct',efficiency_pct, ...
    'VoutEndpointChange_pct',voutChange_pct, ...
    'VaryingPowerMethod',varyingPowerMethod);

harmonicNumber = (1:N_harmonics).';
harmonics = table(harmonicNumber,Hvin.rms,100*Hvin.rms/max(Hvin.rms(1),eps), ...
                  Hiin.rms,100*Hiin.rms/max(Hiin.rms(1),eps), ...
    'VariableNames',{'Harmonic','V_RMS_V','V_pct_of_fundamental', ...
                     'I_RMS_A','I_pct_of_fundamental'});

disp(' ');
disp('First 20 harmonic rows:');
disp(harmonics(1:min(20,height(harmonics)),:));

%% ---------------- OPTIONAL PLOTS ----------------
if makePlots
    figure('Name','Input voltage and current');
    yyaxis left
    plot(tw,vin,'LineWidth',1.0);
    ylabel('Input voltage [V]');
    yyaxis right
    plot(tw,iin,'LineWidth',1.0);
    ylabel('Input current [A]');
    xlabel('Time [s]');
    title(sprintf('Input waveforms: PF = %.4f, current THD = %.2f%%', ...
                  truePF,I_THD_harm_pct));
    grid on;

    figure('Name','Output load powers');
    plot(tw,pFixed,'LineWidth',1.0);
    hold on;
    plot(tw,pVary,'LineWidth',1.0);
    plot(tw,pOut,'LineWidth',1.1);
    hold off;
    xlabel('Time [s]');
    ylabel('Power [W]');
    legend('Fixed load','Varying load','Total useful output', ...
           'Location','best');
    title(sprintf('Average total output power = %.3f W',Pout));
    grid on;

    figure('Name','Input-current harmonics');
    stem(harmonicNumber,100*Hiin.rms/max(Hiin.rms(1),eps),'filled');
    xlabel('Harmonic number');
    ylabel('Current harmonic [% of fundamental]');
    title(sprintf('Input-current harmonics: THD_{2-%d} = %.2f%%', ...
                  N_harmonics,I_THD_harm_pct));
    grid on;
end

%% ================= LOCAL FUNCTIONS =================
function [names,data] = readLtspiceText(fileName)
    fid = fopen(fileName,'r');
    if fid < 0
        error('Cannot open LTspice export: %s',fileName);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    header = fgetl(fid);
    if ~ischar(header)
        error('The LTspice export is empty: %s',fileName);
    end

    if contains(header,sprintf('\t'))
        names = strtrim(strsplit(strtrim(header),sprintf('\t')));
    else
        names = regexp(strtrim(header),'\s+','split');
    end

    nColumns = numel(names);
    formatSpec = repmat('%f',1,nColumns);
    C = textscan(fid,formatSpec, ...
        'Delimiter',sprintf('\t '), ...
        'MultipleDelimsAsOne',true, ...
        'CollectOutput',true, ...
        'ReturnOnError',false);
    data = C{1};

    if isempty(data) || size(data,2) ~= nColumns
        error(['Could not parse the data. The header has %d columns but the ' ...
               'numeric data has %d.'],nColumns,size(data,2));
    end
end

function idx = requireColumn(names,requestedName)
    idx = find(strcmpi(strtrim(names),strtrim(requestedName)),1);
    if isempty(idx)
        fprintf('Available columns:\n');
        for k = 1:numel(names)
            fprintf('  %2d: %s\n',k,names{k});
        end
        error('Required column "%s" was not found.',requestedName);
    end
end

function idx = optionalColumn(names,requestedName)
    idx = find(strcmpi(strtrim(names),strtrim(requestedName)),1);
end

function [tw,Yw] = exactTimeWindow(t,Y,tStart,tStop)
    inside = t > tStart & t < tStop;
    tw = [tStart; t(inside); tStop];
    Ystart = interp1(t,Y,tStart,'linear');
    Ystop  = interp1(t,Y,tStop,'linear');
    Yw = [Ystart; Y(inside,:); Ystop];
end

function value = timeAverage(t,x)
    value = trapz(t,x)/(t(end)-t(1));
end

function value = timeRms(t,x)
    value = sqrt(max(trapz(t,x.^2)/(t(end)-t(1)),0));
end

function H = harmonicAnalysis(t,x,f0,nMax)
    duration = t(end)-t(1);
    H.dc = trapz(t,x)/duration;
    H.rms = zeros(nMax,1);
    H.phase_rad = zeros(nMax,1);
    H.phasorRms = complex(zeros(nMax,1));

    for n = 1:nMax
        omega = 2*pi*n*f0;
        aCos = (2/duration)*trapz(t,x.*cos(omega*t));
        bSin = (2/duration)*trapz(t,x.*sin(omega*t));
        H.phasorRms(n) = (aCos-1i*bSin)/sqrt(2);
        H.rms(n) = abs(H.phasorRms(n));
        H.phase_rad(n) = angle(H.phasorRms(n));
    end
end

function thdRatio = residualTHD(totalRms,dcValue,fundamentalRms)
    residualSquared = totalRms^2-dcValue^2-fundamentalRms^2;
    thdRatio = sqrt(max(residualSquared,0))/max(fundamentalRms,eps);
end
