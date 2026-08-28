clc
clearvars
close all

fs       = 12;
lwLine   = 1;
lwBox    = 3;
fontName = 'Lucida Sans Unicode';

set(groot, 'DefaultAxesFontName',        fontName);
set(groot, 'DefaultAxesFontWeight',      'bold');
set(groot, 'DefaultAxesFontSize',        fs);
set(groot, 'DefaultAxesTitleFontWeight', 'bold');
set(groot, 'DefaultTextFontName',        fontName);
set(groot, 'DefaultTextFontWeight',      'bold');
set(groot, 'DefaultTextFontSize',        fs);
set(groot, 'DefaultLegendFontName',      fontName);
set(groot, 'DefaultLegendFontWeight',    'bold');
set(groot, 'DefaultLegendFontSize',      fs);
set(groot, 'DefaultColorbarFontWeight',  'bold');

st.fs       = fs;
st.lwLine   = lwLine;
st.lwBox    = lwBox;
st.fontName = fontName;
st.mrk      = 8;
st.mrkBig   = 8;

useRounds  = [1 2 3];
showRounds = [4 5];
wlMin      = 840;
wlMax      = 1000;

ana(1).name       = 'Curvature';
ana(1).offset     = 100;
ana(1).dataFolder = 'Data/CURVA977';
ana(1).folderGlob = '*a medicion';
ana(1).fileGlob   = 'sand2cm977nm a';
ana(1).regex      = '^.*CURVA977.(\d+).*medicion.sand2cm977nm a (\d+\.\d+) cm_USB.*$';
ana(1).convFun    = @(v) 2*(v/100) ./ ((v/100).^2 + 0.14^2);
ana(1).selMin     = 0.00;
ana(1).selMax     = 0.51;
ana(1).unit       = 'm^{-1}';
ana(1).appliedLbl = 'Applied curvature [m^{-1}]';
ana(1).estimLbl   = 'Estimated curvature [m^{-1}]';
ana(1).ratioLbl   = 'Curvature [m^{-1}]';
ana(1).matOut     = 'calibracioncurvatura_977.mat';

ana(2).name       = 'Strain';
ana(2).offset     = 200;
ana(2).dataFolder = 'Data/STRAIN977';
ana(2).folderGlob = '*medicion';
ana(2).fileGlob   = 'sand2cm977nm a';
ana(2).regex      = '^.*STRAIN977.(\d+)\D*medicion.sand2cm977nm a (\d+(?:\.\d+)?) um_USB.*$';
ana(2).convFun    = @(v) v/0.04;
ana(2).selMin     = 0;
ana(2).selMax     = Inf;
ana(2).unit       = 'µε';
ana(2).appliedLbl = 'Applied strain [µε]';
ana(2).estimLbl   = 'Estimated strain [µε]';
ana(2).ratioLbl   = 'Strain [µε]';
ana(2).matOut     = 'calibracionstrain_977.mat';

ana(3).name       = 'Temperature';
ana(3).offset     = 300;
ana(3).dataFolder = 'Data/TEMP977';
ana(3).folderGlob = '*a medicion';
ana(3).fileGlob   = 'sand977 a';
ana(3).regex      = '^.*TEMP977.(\d+)\D*a medicion.sand977 a (\d+\.?\d*)\s*.C.*$';
ana(3).convFun    = @(v) v;
ana(3).selMin     = 20;
ana(3).selMax     = 150;
ana(3).unit       = '°C';
ana(3).appliedLbl = 'Applied temperature [°C]';
ana(3).estimLbl   = 'Estimated temperature [°C]';
ana(3).ratioLbl   = 'Temperature [°C]';
ana(3).matOut     = 'calibraciontemperatura_977.mat';

summary = table( ...
    'Size',          [numel(ana) 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Analysis','RMSE','R2','sigma'});

for a = 1:numel(ana)
    r = processAnalysis(ana(a), st, useRounds, showRounds, wlMin, wlMax);

    summary.Analysis(a) = string(ana(a).name);
    summary.RMSE(a)     = r.RMSE;
    summary.R2(a)       = r.R2;
    summary.sigma(a)    = r.sigma;

    fprintf('%-12s | RMSE = %.4f %s | R^2 = %.4f | sigma = %.4f %s\n', ...
        ana(a).name, r.RMSE, ana(a).unit, r.R2, r.sigma, ana(a).unit);

    U       = r.U;       %#ok<NASGU>
    sim     = r.sim;     %#ok<NASGU>
    xpre    = r.xpre;    %#ok<NASGU>
    dsim_dx = r.dsim_dx; %#ok<NASGU>
    wl      = r.wl;      %#ok<NASGU>
    xscale  = r.xscale;  %#ok<NASGU>
    save(ana(a).matOut, 'U','sim','xpre','dsim_dx','wl','xscale');
end

fprintf('\n========== SUMMARY ==========\n');
disp(summary);

set(groot, 'DefaultAxesFontName',        'factory');
set(groot, 'DefaultAxesFontWeight',      'factory');
set(groot, 'DefaultAxesFontSize',        'factory');
set(groot, 'DefaultAxesTitleFontWeight', 'factory');
set(groot, 'DefaultTextFontName',        'factory');
set(groot, 'DefaultTextFontWeight',      'factory');
set(groot, 'DefaultTextFontSize',        'factory');
set(groot, 'DefaultLegendFontName',      'factory');
set(groot, 'DefaultLegendFontWeight',    'factory');
set(groot, 'DefaultLegendFontSize',      'factory');
set(groot, 'DefaultColorbarFontWeight',  'factory');

function r = processAnalysis(cfg, st, useRounds, showRounds, wlMin, wlMax)
    off = cfg.offset;

    folders = dir(fullfile(cfg.dataFolder, cfg.folderGlob));
    if isempty(folders)
        error('No folders matching "%s" in %s', cfg.folderGlob, cfg.dataFolder);
    end

    fdataAll = [];
    wl       = [];
    valAll   = [];
    roundAll = [];

    for k = 1:numel(folders)
        fdir = dir(fullfile(folders(k).folder, folders(k).name, [cfg.fileGlob '*.txt']));
        for j = 1:numel(fdir)
            fname = fullfile(fdir(j).folder, fdir(j).name);
            d = load(fname);
            if isempty(wl)
                wl = d(:,1);
            elseif ~isequal(wl, d(:,1))
                error('Wavelength mismatch in %s', fname);
            end
            tok = regexp(fname, cfg.regex, 'tokens');
            if isempty(tok), continue; end
            fdataAll = [fdataAll, d(:,2)];                       %#ok<AGROW>
            roundAll = [roundAll, str2double(tok{1}{1})];        %#ok<AGROW>
            valAll   = [valAll, cfg.convFun(str2double(tok{1}{2}))]; %#ok<AGROW>
        end
    end

    if isempty(fdataAll)
        error('No files found in %s', cfg.dataFolder);
    end

    valAll   = valAll(:);
    roundAll = roundAll(:);

    sel      = valAll >= cfg.selMin & valAll <= cfg.selMax;
    fdataAll = fdataAll(:, sel);
    valAll   = valAll(sel);
    roundAll = roundAll(sel);

    isCal  = ismember(roundAll, useRounds);
    isShow = ismember(roundAll, showRounds);
    if ~any(isCal)
        error('No calibration rounds present in %s', cfg.name);
    end

    valTest   = valAll(isShow);
    fdataTest = fdataAll(:, isShow);
    mTest     = roundAll(isShow);

    val   = valAll(isCal);
    fdata = fdataAll(:, isCal);
    mCal  = roundAll(isCal);

    idxwl     = wl >= wlMin & wl <= wlMax;
    wl        = wl(idxwl);
    fdata     = fdata(idxwl,:);
    fdataTest = fdataTest(idxwl,:);

    [umed, ~, idxmed] = unique(mCal);
    co = get(groot, 'DefaultAxesColorOrder');

    [U,S,V] = svds(double(fdata), 4);

    ax = newAxes(off+1, st);
    plot(ax, wl, fdata, 'LineWidth', st.lwLine);
    styleAxes(ax, st);
    xlabel(ax, 'Wavelength [nm]'); ylabel(ax, 'Intensity [a.u.]');

    ax = newAxes(off+2, st);
    plot(ax, wl, 10*log10(fdata), 'LineWidth', st.lwLine);
    styleAxes(ax, st);
    xlabel(ax, 'Wavelength [nm]'); ylabel(ax, 'Power [dB]');

    for m = 1:4
        ax = newAxes(off+2+m, st);
        plot(ax, wl, U(:,m), 'LineWidth', st.lwLine, 'Color', co(m,:));
        styleAxes(ax, st);
        xlabel(ax, 'Wavelength [nm]'); ylabel(ax, 'Spectral component');
        legend(ax, sprintf('Mode %d', m), 'Location', 'northwest', 'Box', 'off');
    end

    ampl = V*S;
    ax = newAxes(off+7, st);
    plot(ax, 1:size(ampl,1), ampl, '.-', 'LineWidth', st.lwLine, 'MarkerSize', st.mrk);
    styleAxes(ax, st);
    xlabel(ax, 'Spectrum index'); ylabel(ax, 'Amplitude');
    legend(ax, {'Mode 1','Mode 2','Mode 3','Mode 4'}, 'Location', 'best', 'Box', 'off');
    xlim(ax, [1 size(ampl,1)]);

    ratios = ampl(:,2:4) ./ ampl(:,1);

    axs = newSubplots(off+8, [3 1], st);
    set(gcf, 'Position', [80 40 950 1050]);
    for k = 1:numel(umed)
        s2 = idxmed == k;
        [vs, is] = sort(val(s2));
        a = ratios(s2,:); a = a(is,:);
        for j = 1:3
            plot(axs(j), vs, a(:,j), '.-', 'LineWidth', st.lwLine, ...
                 'MarkerSize', st.mrkBig, 'MarkerEdgeColor', 'k');
        end
    end

    polyorder = 6;
    xscale = max(abs(val));
    xn     = (val/xscale) .^ (0:polyorder);
    cf     = xn \ ratios;

    xpre   = linspace(min(val), max(val), 100)';
    xpre_n = (xpre/xscale) .^ (0:polyorder);
    sim    = xpre_n * cf;

    titlesR = {'A2/A1','A3/A1','A4/A1'};
    for j = 1:3
        plot(axs(j), xpre, sim(:,j), '-', 'LineWidth', st.lwLine);
        title(axs(j), titlesR{j});
        xlabel(axs(j), cfg.ratioLbl);
        set(axs(j), 'XLim', [min(val) max(val)]);
    end
    legend(axs(1), [cellstr(num2str(umed)); {'Fitting'}], 'Location', 'best', 'Box', 'off');

    dsim_dx = (xpre_n(:,1:end-1) * ((1:polyorder)' .* cf(2:end,:))) / xscale;

    figure(off+9); clf
    set(gcf, 'Position', [100 80 950 850]);
    ax5 = axes('Parent', gcf, 'NextPlot', 'add');
    set(ax5, 'FontSize', st.fs, 'FontWeight', 'bold', 'FontName', st.fontName, ...
             'LineWidth', st.lwBox, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    box(ax5, 'on'); grid(ax5, 'on'); view(ax5, 3);
    hExp = gobjects(numel(umed),1);
    for k = 1:numel(umed)
        s2 = idxmed == k;
        [~, is] = sort(val(s2));
        a = ratios(s2,:); a = a(is,:);
        hExp(k) = plot3(ax5, a(:,1), a(:,2), a(:,3), '.-', 'LineWidth', st.lwLine, ...
                        'MarkerSize', st.mrkBig, 'MarkerEdgeColor', 'k');
    end
    hFit = plot3(ax5, sim(:,1), sim(:,2), sim(:,3), 'k-', 'LineWidth', st.lwLine);
    xlabel(ax5, 'A2/A1'); ylabel(ax5, 'A3/A1'); zlabel(ax5, 'A4/A1');
    legend(ax5, [hExp(1); hFit], {'Experimental data Tests','Fitting'}, ...
           'Location', 'best', 'Box', 'off');

    xEstCal = demodulate((U'*fdata)', ratios, sim, xpre, dsim_dx);

    if ~isempty(valTest)
        xEstTest = demodulate((U'*fdataTest)', ratios, sim, xpre, dsim_dx);
    end

    p     = polyfit(val, xEstCal, 1);
    resid = xEstCal - polyval(p, val);
    R2    = 1 - sum(resid.^2) / sum((xEstCal - mean(xEstCal)).^2);
    sigma = std(resid);

    if ~isempty(valTest)
        RMSE = sqrt(mean((xEstTest - valTest).^2));
    else
        RMSE = sqrt(mean((xEstCal - val).^2));
    end

    if ~isempty(valTest)
        vPlot = valTest; ePlot = xEstTest;
        [umedP, ~, idxmedP] = unique(mTest);
        pref = 'Test ';
    else
        vPlot = val; ePlot = xEstCal;
        umedP = umed; idxmedP = idxmed;
        pref = 'Round ';
    end

    ax = newAxes(off+10, st);
    h    = gobjects(numel(umedP),1);
    labs = cell(numel(umedP),1);
    for k = 1:numel(umedP)
        s2 = idxmedP == k;
        [vs, is] = sort(vPlot(s2));
        es = ePlot(s2); es = es(is);
        h(k) = plot(ax, vs, es, '.-', 'Color', co(mod(k-1,size(co,1))+1,:), ...
             'LineWidth', st.lwLine, 'MarkerSize', st.mrkBig, 'MarkerEdgeColor', 'k');
        labs{k} = sprintf('%s%d', pref, umedP(k));
    end
    xln = linspace(min(vPlot), max(vPlot), 200);
    hf  = plot(ax, xln, polyval(p, xln), 'k--', 'LineWidth', 1.5);
    xlim(ax, [min(vPlot) max(vPlot)]);
    ylim(ax, [min(ePlot) max(ePlot)]);
    text(ax, 0.05, 0.90, ...
        {sprintf('R^2 = %.4f', R2), sprintf('\\sigma = %.4f %s', sigma, cfg.unit)}, ...
        'Units', 'normalized', 'BackgroundColor', 'white', 'EdgeColor', 'none');
    styleAxes(ax, st);
    xlabel(ax, cfg.appliedLbl); ylabel(ax, cfg.estimLbl);
    legend(ax, [h; hf], [labs; {'Linear fit'}], 'Location', 'southeast', 'Box', 'off');

    r.RMSE    = RMSE;
    r.R2      = R2;
    r.sigma   = sigma;
    r.U       = U;
    r.sim     = sim;
    r.xpre    = xpre;
    r.dsim_dx = dsim_dx;
    r.wl      = wl;
    r.xscale  = xscale;
end

function xest = demodulate(amplNew, ratiosCal, sim, xpre, dsim_dx)  %#ok<INUSD>
    rat  = amplNew(:,2:4) ./ amplNew(:,1);
    dist = squeeze(sum((rat - permute(sim, [3 2 1])).^2, 2));
    [~, imin] = min(dist, [], 2);
    xrgh  = xpre(imin);
    delta = rat - sim(imin,:);
    xcorr = sum(dsim_dx(imin,:) .* delta, 2) ./ sum(dsim_dx(imin,:).^2, 2);
    xest  = xrgh + xcorr;
end

function ax = newAxes(fignum, st)
    figure(fignum); clf
    ax = gca;
    set(ax, 'FontSize', st.fs, 'FontWeight', 'bold', 'FontName', st.fontName);
    hold(ax, 'on');
end

function styleAxes(ax, st)
    ax.LineWidth = st.lwBox;
    ax.XColor = 'k';
    ax.YColor = 'k';
    try, ax.ZColor = 'k'; catch, end
    ax.FontSize = st.fs;
    ax.FontWeight = 'bold';
    box(ax, 'on'); grid(ax, 'off');
end

function axs = newSubplots(fignum, layout, st)
    figure(fignum); clf
    n = prod(layout);
    axs = gobjects(n,1);
    for k = 1:n
        axs(k) = subplot(layout(1), layout(2), k);
        set(axs(k), 'FontSize', st.fs, 'FontWeight', 'bold', 'FontName', st.fontName, ...
                    'LineWidth', st.lwBox, 'XColor', 'k', 'YColor', 'k');
        box(axs(k), 'on'); grid(axs(k), 'off'); hold(axs(k), 'on');
    end
end