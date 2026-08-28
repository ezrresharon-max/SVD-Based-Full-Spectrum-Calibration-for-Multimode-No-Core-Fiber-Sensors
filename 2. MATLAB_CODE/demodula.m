N_COMP = 3;

ROLES = {'Curvature','Strain','Temperature'};
UNIT  = {'m^-1',      'µε',    '°C'};

CAND_U      = {'U'};
CAND_SIM    = {'sim','Sim','S_sim'};
CAND_PARAM  = {'xpre','cmpre','Tpre','strain','eps','T','temp','param','x','vals'};
CAND_DPARAM = {'dsim_dx','dsimdx','dSdx','jacobian','J'};
CAND_WL     = {'wl','wavelength','lambda','lam'};

TRUE_REGEX = {'(\d+\.?\d*)\s*cm_USB', '(\d+\.?\d*)\s*um_USB', '(\d+\.?\d*)\s*.C'};
TRUE_CONV  = {@(v) 2.*(v/100)./((v/100).^2 + 0.14^2), @(v) v/0.04, @(v) v};

nR = numel(ROLES);

for r = 1:nR
    [f,p] = uigetfile('*.mat', sprintf('Calibration for %s', ROLES{r}));
    if isequal(f,0), error('No calibration selected for %s.', ROLES{r}); end
    d = load(fullfile(p,f));

    [~, U]        = find_first(d, CAND_U);
    [~, sim]      = find_first(d, CAND_SIM);
    [nP, axisVals]= find_first(d, CAND_PARAM);
    [~, dsim_dx]  = find_first(d, CAND_DPARAM);
    [~, wl]       = find_first(d, CAND_WL);

    miss = {};
    if isempty(U),        miss{end+1} = 'U'; end %#ok<*SAGROW>
    if isempty(sim),      miss{end+1} = 'sim'; end
    if isempty(axisVals), miss{end+1} = 'axis (xpre/cmpre/Tpre)'; end
    if isempty(dsim_dx),  miss{end+1} = 'dsim_dx'; end
    if isempty(wl),       miss{end+1} = 'wl'; end
    if ~isempty(miss)
        error('Calibration %s: missing %s. Present: %s', ...
            ROLES{r}, strjoin(miss,', '), strjoin(fieldnames(d),', '));
    end

    cal(r).label    = ROLES{r};
    cal(r).unit     = UNIT{r};
    cal(r).U        = U;
    cal(r).sim      = sim;
    cal(r).axisVals = axisVals(:);
    cal(r).dsim_dx  = dsim_dx;
    cal(r).wl       = wl(:);
    cal(r).Nc       = min([N_COMP, size(sim,2), size(dsim_dx,2)]);

    S = sim(:,1:cal(r).Nc);
    cal(r).radius = sqrt(mean(sum((S - mean(S,1)).^2, 2)));
    if cal(r).radius < eps, cal(r).radius = 1; end

    fprintf('OK  %-12s <- %-32s  [axis=%s, Nc=%d, range=%.4g..%.4g %s]\n', ...
        cal(r).label, f, nP, cal(r).Nc, min(cal(r).axisVals), max(cal(r).axisVals), cal(r).unit);
end
fprintf('\n');

[files, pathspec] = uigetfile( ...
    {'*.txt;*.dat;*.csv','Spectra'}, ...
    'Select spectrum(a) to identify', p, 'MultiSelect','on');
if isequal(files,0), error('No spectrum selected.'); end
if ischar(files), files = {files}; end
Nsp = numel(files);

Winner    = strings(Nsp,1);
EstValue  = nan(Nsp,1);
WinUnit   = strings(Nsp,1);
TrueValue = nan(Nsp,1);
ErrValue  = nan(Nsp,1);

errStore = cell(nR,1);

for s = 1:Nsp
    name = files{s};

    try
        raw = load(fullfile(pathspec,name));
    catch
        tmp = importdata(fullfile(pathspec,name));
        if isstruct(tmp), raw = tmp.data; else, raw = tmp; end
    end
    wl_raw = raw(:,1); spec_raw = raw(:,2);
    ok = isfinite(wl_raw) & isfinite(spec_raw);
    wl_raw = wl_raw(ok); spec_raw = spec_raw(ok);
    [wl_raw, ord] = sort(wl_raw); spec_raw = spec_raw(ord);

    for r = 1:nR
        Nc  = cal(r).Nc;
        wlk = cal(r).wl;

        idx = wl_raw >= min(wlk) & wl_raw <= max(wlk);
        if nnz(idx) < 10
            cal(r).resNorm = Inf; cal(r).value = NaN;
            cal(r).idxmin = NaN; cal(r).ampln = nan(Nc,1); continue
        end
        spec_i = interp1(wl_raw(idx), spec_raw(idx), wlk, 'linear', 'extrap');

        ampl = cal(r).U' * spec_i;
        if abs(ampl(1)) < eps
            cal(r).resNorm = Inf; cal(r).value = NaN;
            cal(r).idxmin = NaN; cal(r).ampln = nan(Nc,1); continue
        end
        ampln = ampl(2:1+Nc) ./ ampl(1);

        dif  = cal(r).sim(:,1:Nc) - ampln.';
        dist = sum(dif.^2, 2);
        [distMin, idxmin] = min(dist);

        delta = ampln.' - cal(r).sim(idxmin,1:Nc);
        grad  = cal(r).dsim_dx(idxmin,1:Nc);
        denom = sum(grad.^2);
        if denom > 0, corr = sum(grad.*delta)/denom; else, corr = 0; end

        cal(r).value   = cal(r).axisVals(idxmin) + corr;
        cal(r).resNorm = sqrt(distMin) / cal(r).radius;
        cal(r).idxmin  = idxmin;
        cal(r).ampln   = ampln;
    end

    rn = [cal.resNorm];
    [~, win] = min(rn);

    trueRole = 0; trueVal = NaN;
    for r = 1:nR
        tk = regexp(name, TRUE_REGEX{r}, 'tokens', 'once');
        if ~isempty(tk)
            trueRole = r; trueVal = TRUE_CONV{r}(str2double(tk{1})); break
        end
    end

    Winner(s)    = cal(win).label;
    EstValue(s)  = cal(win).value;
    WinUnit(s)   = cal(win).unit;
    TrueValue(s) = trueVal;

    fprintf('%s\n', name);
    if trueRole == 0
        fprintf('  Winner: %s = %.4g %s\n', cal(win).label, cal(win).value, cal(win).unit);
    elseif trueRole == win
        err = cal(win).value - trueVal;
        ErrValue(s) = err;
        errStore{win} = [errStore{win}; err];
        fprintf('  Winner: %s = %.4g %s  | true = %.4g %s  | error = %.4g %s\n', ...
            cal(win).label, cal(win).value, cal(win).unit, trueVal, cal(win).unit, err, cal(win).unit);
    else
        fprintf('  Winner: %s = %.4g %s  | MISCLASSIFIED (true: %s = %.4g %s)\n', ...
            cal(win).label, cal(win).value, cal(win).unit, ...
            cal(trueRole).label, trueVal, cal(trueRole).unit);
    end
end

fprintf('\n========== ERROR VS TRUE VALUE ==========\n');
for r = 1:nR
    e = errStore{r};
    if isempty(e), continue; end
    fprintf('%-12s N=%d  RMSE=%.4f %s  MAE=%.4f %s  std=%.4f %s\n', ...
        cal(r).label, numel(e), sqrt(mean(e.^2)), cal(r).unit, ...
        mean(abs(e)), cal(r).unit, std(e), cal(r).unit);
end

Summary = table(string(files(:)), Winner, EstValue, WinUnit, TrueValue, ErrValue, ...
    'VariableNames', {'File','Parameter','Value','Unit','TrueValue','Error'});
disp(Summary)
assignin('base','Summary',Summary);

colors = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19];

figure('Color','w','Position',[120 200 1050 400]);
subplot(1,2,1)
rn = [cal.resNorm];
b = bar(1:nR, rn); b.FaceColor = 'flat';
for r = 1:nR
    if r==win, b.CData(r,:) = [0.20 0.70 0.20]; else, b.CData(r,:) = [0.78 0.78 0.78]; end
end
xticks(1:nR); xticklabels({cal.label}); ylabel('Normalized residual (lower = better)')
title('Fit per calibration'); grid on

subplot(1,2,2)
for r = 1:nR
    if isfinite(cal(r).resNorm) && cal(r).Nc >= 2
        plot(cal(r).sim(:,1), cal(r).sim(:,2), '.-', 'Color', colors(r,:), ...
             'MarkerSize', 7, 'DisplayName', cal(r).label); hold on
        plot(cal(r).ampln(1), cal(r).ampln(2), 'x', 'Color', colors(r,:), ...
             'MarkerSize', 13, 'LineWidth', 2, 'HandleVisibility','off')
    end
end
xlabel('A2/A1'); ylabel('A3/A1'); grid on; axis equal
title('Manifolds and measured point (x)'); legend('Location','best')

function [name, value] = find_first(d, candidates)
    name = ''; value = [];
    for j = 1:numel(candidates)
        if isfield(d, candidates{j})
            name  = candidates{j};
            value = d.(candidates{j});
            return
        end
    end
end