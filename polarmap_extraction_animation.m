%% rb_polarmap_anim.m
% Animated simulation of the effect of 82Rb extraction and extraction
% correction on myocardial blood flow (MBF) polar maps.
%
% The MBF in a mid-anterior defect is swept from healthy (MBFh) down to 0.5,
% held, then brought back up to MBFd_end. Five maps are shown side by side:
% noisy [15O]H2O MBF, 82Rb retention, 82Rb K1, and the MBF recovered from
% retention and from K1 by inverting the respective extraction functions.
% A sixth panel shows the defect/healthy contrast for each tracer.
%
% Output:  GIF + MP4
%
% Requirements
%   - MATLAB with Image Processing Toolbox (imgaussfilt) and Statistics
%     Toolbox (prctile)
%   - MATLAB R2019b or newer (tiledlayout, local functions in scripts)
% 
% Author:  Lars Tolbod, Aarhus University, 2026

clear; close all

%% ---------- Settings ----------
MBFh     = 2.5;     % MBF in healthy myocardium [mL/min/g]
MBFd_end = 1.8;     % MBF in the defect at the end of the animation [mL/min/g]

% Defect MBF per frame: healthy -> 0.5 (90 frames), hold (25), -> MBFd_end (50), hold (25)
Fd_seq = [linspace(MBFh, 0.5, 90),      0.5*ones(1, 25), ...
          linspace(0.5, MBFd_end, 50),  MBFd_end*ones(1, 25)];

cv       = 0.10;    % relative noise (CoV) on each measured map: H2O MBF, Rb retention, Rb K1
smoothPx = 1.5;     % spatial correlation of the noise and blurring of the true flow [pixels]
showContrast = true;    % plot defect contrast in the empty bottom-left tile
fps      = 15;
outGif   = 'rb_polarmap_anim.gif';
outMp4   = 'rb_polarmap_anim.mp4';

% Extraction models (parameters as in the static analysis)
pRet = [0.74 0.083];                                % retention: R = (1 - exp(-(a + b*F)/F)) * F
pK1  = [0.77 0.63];                                 % K1:        K1 = (1 - a*exp(-b/F)) * F
fRet = @(F) (1 - exp(-(pRet(1) + pRet(2).*F)./F)).*F;
fK1  = @(F) (1 - pK1(1)*exp(-pK1(2)./F)).*F;

% The retention function has no closed-form inverse -> lookup table
Fgrid   = linspace(0.01, 15, 20000);
retGrid = fRet(Fgrid);                              % monotonically increasing
invRet  = @(R) interp1(retGrid, Fgrid, R, 'linear', 'extrap');

%% ---------- Input data ----------
% Polar-map pixel grid: 28 rings (radius 28 at row 1 down to 1 at row 28) by
% 73 angular columns (0:5:360 degrees, closing column included).
% Convention: X = -r*cos(theta), Y = r*sin(theta), consistent with the
% territory lines drawn in plotPolarmapLines.
sz    = [28 73];                                    % polar-map size [radial x angular]
r     = (sz(1):-1:1)';                              % ring radii [pixels]
theta = 0:5:360;                                    % angles [degrees]
data.rest.polarmap.X = -r * cosd(theta);
data.rest.polarmap.Y =  r * sind(theta);

% Coronary-territory line settings
data.rest.TAC.exclude      = 1;                     % border exclusion [10-degree steps]  (SET TO YOUR VALUE)
data.rest.TAC.exclude_apex = true;                  % exclude apex ring from territory lines (SET TO YOUR VALUE)

%% ---------- Segment index grid (17-segment-style model) ----------
nAng = 72;                                          % angular samples
nRad = 28;                                          % radial samples
indices = zeros(nAng, nRad);

basalCols = 1 : 2*nRad/7;                           % basal ring
midCols   = 2*nRad/7+1 : 4*nRad/7;                  % mid ring
apexCols  = 4*nRad/7+1 : 6*nRad/7;                  % apical ring

% Six 60-degree sectors per ring; segment numbers in sector order
% (antero-septal, anterior, antero-lateral, infero-lateral, inferior, infero-septal)
basalSeg = [6 1 2 3 4 5];
midSeg   = [12 7 8 9 10 11];
for s = 1:6
    rows = (s-1)*nAng/6+1 : s*nAng/6;
    indices(rows, basalCols) = basalSeg(s);
    indices(rows, midCols)   = midSeg(s);
end

% Apex: three 90-degree sectors (anterior, lateral, inferior), rotated by nAng/8
for s = 1:3
    rows = (s*nAng/4+1 : (s+1)*nAng/4) - round(nAng/8);
    indices(rows, apexCols) = 12 + s;               % segments 13, 14, 15
end

% Defect = mid anterior (segment 7). The polar map has one more angular
% column than the index grid (closing column), so pad the mask to full size.
m = indices' == 7;
defectMask = false(sz);
defectMask(1:size(m,1), 1:size(m,2)) = m;

%% ---------- Fixed noise realisation (no flicker between frames) ----------
% Gaussian noise, smoothed for spatial correlation, then rescaled to zero
% mean / unit SD so the CoV in the displayed measured maps equals cv.
rng(1)
U = cell(1, 3);
for k = 1:3
    z    = imgaussfilt(randn(sz), smoothPx);
    U{k} = (z - mean(z(:))) ./ std(z(:));
end
healthy = imgaussfilt(double(defectMask), 2*smoothPx) < 0.01;   % pixels well away from the defect

%% ---------- Figure ----------
% Layout:  row 1: H2O MBF  | Rb retention         | Rb K1
%          row 2: contrast | Rb MBF (retention)   | Rb MBF (K1)
names = {'[^{15}O]H_2O MBF', '^{82}Rb Retention', '^{82}Rb K_1', ...
         '^{82}Rb MBF (from retention)', '^{82}Rb MBF (from K_1)'};
nMap    = numel(names);
tilePos = [1 2 3 5 6];

fig = figure('Position', [10 10 1300 820], 'Color', 'w');
tl  = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Colour limits are fixed from the no-defect maps so the scale does not jump
M0   = simMaps(MBFh, MBFh, defectMask, U, cv, smoothPx, fRet, fK1, pK1, invRet);
cmax = zeros(1, nMap);
for k = 1:nMap
    cmax(k) = prctile(M0{k}(:), 90);
end
cmax([1 4 5]) = max(cmax([1 4 5]));                 % common scale on the three MBF maps

M  = simMaps(Fd_seq(1), MBFh, defectMask, U, cv, smoothPx, fRet, fK1, pK1, invRet);
ax = cell(1, nMap);
h  = cell(1, nMap);
for k = 1:nMap
    ax{k} = nexttile(tl, tilePos(k));
    h{k}  = pcolor(ax{k}, data.rest.polarmap.X, data.rest.polarmap.Y, M{k});
    colorbar(ax{k}, 'Color', 'b');
    shading(ax{k}, 'flat');
    set(ax{k}, 'DataAspectRatio', [1 1 1]);
    caxis(ax{k}, [0 cmax(k)]);                      
    axis(ax{k}, 'tight');
    axis(ax{k}, 'off');
    colormap(ax{k}, get_germano);
    plotPolarmapLines(ax{k}, data);
    title(ax{k}, mapTitle(names{k}, M{k}, healthy));
end
sg = sgtitle(fig, '');

% Contrast panel: defect/healthy ratio for each map (noise-free curves)
if showContrast
    axC     = nexttile(tl, 4);
    Fd_axis = linspace(0.3, MBFh, 300);
    rC  = {Fd_axis/MBFh, fRet(Fd_axis)/fRet(MBFh), fK1(Fd_axis)/fK1(MBFh)};
    col = {[0.3 0.3 0.3], [0.82 0.29 0.36], [0.12 0.37 0.62]};
    hold(axC, 'on')
    hl = gobjects(1, 3);
    hd = gobjects(1, 3);
    for k = 1:3
        hl(k) = plot(axC, Fd_axis, rC{k}, 'Color', col{k}, 'LineWidth', 2.5);
        hd(k) = plot(axC, Fd_seq(1), 1, 'o', 'Color', col{k}, ...
                     'MarkerFaceColor', col{k}, 'MarkerSize', 8);
    end
    xl = xline(axC, Fd_seq(1), ':k');
    hold(axC, 'off')
    set(axC, 'XDir', 'reverse', 'Box', 'off', 'FontSize', 10);
    xlim(axC, [0.3 MBFh]);
    ylim(axC, [0 1.05]);
    xlabel(axC, 'Defect MBF [mL/min/g]');
    ylabel(axC, 'Defect / healthy');
    title(axC, 'Defect contrast')
    legend(axC, hl, {'[^{15}O]H_2O', '^{82}Rb K_1', '^{82}Rb Retention'}, ...
           'Location', 'southwest', 'Box', 'off', 'FontSize', 9);
    txt = text(axC, 0.97, 0.97, '', 'Units', 'normalized', ...
               'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);
end

%% ---------- Animate ----------
v = VideoWriter(outMp4, 'MPEG-4');
v.FrameRate = fps;
v.Quality   = 95;
open(v);

for f = 1:numel(Fd_seq)
    Fd = Fd_seq(f);
    M  = simMaps(Fd, MBFh, defectMask, U, cv, smoothPx, fRet, fK1, pK1, invRet);

    for k = 1:nMap
        set(h{k}, 'CData', M{k});
        title(ax{k}, mapTitle(names{k}, M{k}, healthy));
    end
    sg.String = sprintf('MBF_{healthy} = %.2f, MBF_{defect} = %.2f', MBFh, Fd);

    if showContrast
        r = [Fd/MBFh, fRet(Fd)/fRet(MBFh), fK1(Fd)/fK1(MBFh)];   % H2O, Rb Ret, Rb K1
        for k = 1:3
            set(hd(k), 'XData', Fd, 'YData', r(k));
        end
        xl.Value   = Fd;
        txt.String = {sprintf('True (H_2O): %.1f : 1',   1/r(1)), ...
                      sprintf('^{82}Rb K_1: %.1f : 1',   1/r(3)), ...
                      sprintf('^{82}Rb Ret: %.1f : 1',   1/r(2))};
    end

    drawnow
    fr  = getframe(fig);
    img = fr.cdata(1:2*floor(end/2), 1:2*floor(end/2), :);      % even size required by MPEG-4
    writeVideo(v, img);

    [A, map] = rgb2ind(img, 256);
    if f == 1
        imwrite(A, map, outGif, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
    else
        imwrite(A, map, outGif, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end
close(v);

%% ---------- Local functions ----------
function M = simMaps(Fd, MBFh, defectMask, U, cv, smoothPx, fRet, fK1, pK1, invRet)
% Simulate the five displayed maps for a given defect MBF Fd.
    F = MBFh * ones(size(defectMask));
    F(defectMask) = Fd;
    F = imgaussfilt(F, smoothPx);                   % finite resolution: soft defect edges

    mbfN = F      .* (1 + cv.*U{1});                % same relative noise on every measurement
    rEqN = fRet(F).* (1 + cv.*U{2});
    K1N  = fK1(F) .* (1 + cv.*U{3});

    mbfRet = invRet(rEqN);                          % MBF recovered from retention
    mbfK1  = lortiecalc(K1N, pK1);                  % MBF recovered from K1
    M = {mbfN, rEqN, K1N, mbfRet, mbfK1};
end

function s = mapTitle(name, map, healthy)
% Panel title with the CoV measured in the healthy myocardium.
    s = sprintf('%s  (CoV %.0f%%)', name, 100*std(map(healthy))/mean(map(healthy)));
end

function plotPolarmapLines(ax, data)
% Draw ring/segment borders (black) and coronary territories on a polar map:
% LAD (red), LCx (blue), RCA (green).
    hold(ax, 'on')
    nR = size(data.rest.polarmap.X, 1);
    r1 = nR/7;  r2 = 3*nR/7;  r3 = 5*nR/7;  r4 = nR;
    a  = 0:360;

    % Rings
    for r = [r1 r2 r3 r4]
        plot(ax, r*cosd(a), r*sind(a), 'k')
    end
    % Radial borders: basal/mid (60-degree spacing) and apical (90-degree spacing)
    for i = 0:60:300
        plot(ax, [r2 r4]*cosd(i), [r2 r4]*sind(i), 'k')
    end
    for i = 45:90:360
        plot(ax, [r2 r1]*cosd(i), [r2 r1]*sind(i), 'k')
    end

    ex = data.rest.TAC.exclude;                     % border exclusion, in 10-degree steps
    e  = ex*10;                                     % border exclusion in degrees

    %% LAD (red)
    L = [];
    L = [L; -r4*cosd(ex*(10:10:(120-e))'), r4*sind(ex*(10:10:(120-e))')];
    L = [L; -[r4 r2]'*cosd(120-e),         [r4 r2]'*sind(120-e)];
    L = [L; -r2*cosd((120:1:135)'-e),      r2*sind((120:1:135)'-e)];
    L = [L; -[r2 r1]'*cosd(135-e),         [r2 r1]'*sind(135-e)];
    if data.rest.TAC.exclude_apex
        L = [L; -r1*cosd((135-e:-1:-45+e)'), r1*sind((135-e:-1:-45+e)')];
    else
        L = [L; -r1*cosd((135-e:1:315+e)'),  r1*sind((135-e:1:315+e)')];
    end
    L = [L; -[r1 r2]'*cosd(-45+e),         [r1 r2]'*sind(-45+e)];
    L = [L; -r2*cosd((-45:1:0)'+e),        r2*sind((-45:1:0)'+e)];
    L = [L; -[r2 r4]'*cosd(e),             [r2 r4]'*sind(e)];
    plot(ax, L(:,1), L(:,2), 'r', 'LineWidth', 2)

    %% LCx (blue)
    L = [];
    L = [L; -r4*cosd(120+(e:10:(120-e))'),     r4*sind(120+(e:10:(120-e))')];
    L = [L; -[r4 r2]'*cosd(240-e),             [r4 r2]'*sind(240-e)];
    L = [L; -r2*cosd((225:1:240)'-e),          r2*sind((225:1:240)'-e)];
    L = [L; -[r2 r1]'*cosd(225-e),             [r2 r1]'*sind(225-e)];
    L = [L; -r1*cosd(135+((90-e):-10:e)'),     r1*sind(135+((90-e):-10:e)')];
    L = [L; -[r1 r2]'*cosd(135+e),             [r1 r2]'*sind(135+e)];
    L = [L; -r2*cosd((135:-1:120)'+e),         r2*sind((135:-1:120)'+e)];
    L = [L; -[r2 r4]'*cosd(120+e),             [r2 r4]'*sind(120+e)];
    plot(ax, L(:,1), L(:,2), 'b', 'LineWidth', 2)

    %% RCA (green)
    L = [];
    L = [L; -r4*cosd(240+(e:10:(120-e))'),     r4*sind(240+(e:10:(120-e))')];
    L = [L; -[r4 r2]'*cosd(360-e),             [r4 r2]'*sind(360-e)];
    L = [L; -r2*cosd((360:-1:315)'-e),         r2*sind((360:-1:315)'-e)];
    L = [L; -[r2 r1]'*cosd(315-e),             [r2 r1]'*sind(315-e)];
    L = [L; -r1*cosd(225+((90-e):-10:e)'),     r1*sind(225+((90-e):-10:e)')];
    L = [L; -[r1 r2]'*cosd(225+e),             [r1 r2]'*sind(225+e)];
    L = [L; -r2*cosd((225:1:240)'+e),          r2*sind((225:1:240)'+e)];
    L = [L; -[r2 r4]'*cosd(240+e),             [r2 r4]'*sind(240+e)];
    plot(ax, L(:,1), L(:,2), 'g', 'LineWidth', 2)

    hold(ax, 'off')
end

%% ---------- Local function: MBF from 82Rb K1 (inverse extraction function) ----------
function F = lortiecalc(K1, p)
%LORTIECALC  Invert the 82Rb K1 extraction function to obtain MBF.
%
%   F = LORTIECALC(K1, p) returns the myocardial blood flow F [mL/min/g]
%   that corresponds to the 82Rb uptake rate constant K1 [mL/min/g],
%   assuming the extraction model
%
%       K1 = (1 - p(1)*exp(-p(2)/F)) * F
%
%   (Lortie et al., Eur J Nucl Med Mol Imaging 2007). K1 may be a scalar,
%   vector or matrix (e.g. a whole polar map); F has the same size.
%
%   Example (parameters used in the Rb polar-map animation):
%       F = lortiecalc(K1map, [0.77 0.63]);
%
%   Author: Lars Tolbod, Aarhus University, 2026

    fK1    = @(F) (1 - p(1)*exp(-p(2)./F)).*F;
    Fgrid  = linspace(0.01, 20, 20000);     % step 0.001 mL/min/g
    K1grid = fK1(Fgrid);                    % monotonically increasing

    F = interp1(K1grid, Fgrid, K1, 'linear', 'extrap');
end

%% ---------- Local function: polar-map colour map ----------
function cmap = get_germano
%GET_GERMANO  64-level colour map used for the polar maps.
%
%   cmap = GET_GERMANO returns a 64-by-3 RGB matrix (values 0-1) that can be
%   passed directly to COLORMAP, e.g.
%
%       colormap(ax, get_germano)
%
%   The scale runs black -> teal -> purple -> orange -> pale yellow-white.

    cmap = [ ...
             0         0         0
        0.0002    0.0309    0.0262
        0.0005    0.0618    0.0525
        0.0007    0.0926    0.0787
        0.0010    0.1235    0.1049
        0.0012    0.1544    0.1311
        0.0015    0.1853    0.1574
        0.0017    0.2162    0.1836
        0.0020    0.2471    0.2098
        0.0022    0.2779    0.2360
        0.0025    0.3088    0.2623
        0.0027    0.3397    0.2885
        0.0029    0.3706    0.3147
        0.0032    0.4015    0.3409
        0.0034    0.4324    0.3672
        0.0037    0.4632    0.3934
        0.0039    0.4941    0.4196
        0.0332    0.4651    0.4245
        0.0625    0.4360    0.4293
        0.0918    0.4069    0.4341
        0.1211    0.3779    0.4390
        0.1504    0.3488    0.4438
        0.1797    0.3197    0.4487
        0.2090    0.2907    0.4535
        0.2383    0.2616    0.4584
        0.2676    0.2325    0.4632
        0.2969    0.2035    0.4681
        0.3262    0.1744    0.4729
        0.3555    0.1453    0.4777
        0.3848    0.1163    0.4826
        0.4141    0.0872    0.4874
        0.4434    0.0581    0.4923
        0.4727    0.0291    0.4971
        0.5020         0    0.5020
        0.5121    0.0135    0.4908
        0.5223    0.0269    0.4797
        0.5325    0.0404    0.4685
        0.5427    0.0539    0.4573
        0.5529    0.0674    0.4462
        0.5842    0.1088    0.4119
        0.6155    0.1503    0.3775
        0.6469    0.1917    0.3432
        0.6782    0.2332    0.3089
        0.7095    0.2746    0.2746
        0.7409    0.3161    0.2403
        0.7722    0.3575    0.2059
        0.8035    0.3990    0.1716
        0.8348    0.4404    0.1373
        0.8662    0.4819    0.1030
        0.8975    0.5233    0.0686
        0.9288    0.5648    0.0343
        0.9602    0.6063         0
        0.9617    0.6386    0.0829
        0.9632    0.6710    0.1658
        0.9647    0.7034    0.2487
        0.9661    0.7357    0.3316
        0.9676    0.7681    0.4145
        0.9691    0.8005    0.4973
        0.9706    0.8328    0.5802
        0.9721    0.8652    0.6631
        0.9736    0.8976    0.7460
        0.9751    0.9299    0.8289
        0.9766    0.9623    0.9118
        0.9781    0.9947    0.9947];
end