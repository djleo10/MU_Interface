function plotSingleMUAP(MUAP, MU_ID, Fsamp, BackColor, FrontColor, crossCorrMatrix)
    Fontssize1 = 10;
    YLimMin = 10000;
    YLimMax = 0;
    col{1} = FrontColor; %[0.5,0.5,1.0];
    col{2} = FrontColor; % [1,1,0];
    YLimCorrect = 0.9;

    h = figure;
    set(h, 'PaperPositionMode', 'auto');
    Position = get(0, 'ScreenSize');
    set(h, 'Position', Position);
    set(h, 'Color', BackColor);

    rY = size(MUAP, 1);
    cY = size(MUAP, 2);

    for r = 1:rY
        for c = 1:cY
            if ~isempty(MUAP{r, c})
                T = floor(size(MUAP{r, c}, 2) / 2);
            end
        end
    end

    tmp = [];
    for r = 1:rY
        for c = 1:cY
            mMUAPs = [];
            if isempty(MUAP{r, c})
                tmp(end + 1) = 0;
            else
                for k1 = 1:size(MUAP{r, c}, 1)
                    [mMUAPs(k1), iMUAPs] = max(abs(MUAP{r, c}(k1, :)));
                end
                tmp(end + 1) = 1.0 * max(mMUAPs) + eps;
            end
            tmp(end) = tmp(end) + eps;
        end
    end

    for r = 1:rY
        dh(r) = 1.5 * max(tmp);
    end

    yLineHight = round(max(tmp) * 100000) / 100;
    if yLineHight > 10
        yLineHight = floor(max(tmp) / 10) * 10;
    end
    xLineWidth = round(0.02 * Fsamp);

    for r = 1:rY
        for c = 1:cY
            if ~isempty(MUAP{r, c})
                mMUAPs = [];
                if size(MUAP{r, c}, 1) > 1
                    for k1 = 1:size(MUAP{r, c}, 1)
                        plot((c - 1) * T * 1.1 + [1:T], MUAP{r, c}(k1, :)' + sum(dh(1:r)), 'Color', [0.6, 0.6, 0.6]);
                    end
                end
                plot((c - 1) * (2 * T + 1) * 1.1 + [1:(2 * T + 1)], mean(MUAP{r, c}, 1)' + sum(dh(1:r)), 'LineWidth', 0.5, 'Color', [0.1, 0.1, 0.7]);
                [mMUAPs(k1), iMUAPs] = max(abs(mean(MUAP{r, c}, 1)));
            end

            YLimMax = max([YLimMax, max(mMUAPs) + sum(dh(1:r))]);
            YLimMin = min([YLimMin, min(mMUAPs) + sum(dh(1:r))]);

            % Add cross-correlation scores
            if r < rY && ~isempty(crossCorrMatrix)
                text((c - 1) * (2 * T + 1) * 1.1 + T, sum(dh(1:r)) + max(mMUAPs) / 2, num2str(crossCorrMatrix(r, c)), 'Color', 'red', 'FontSize', Fontssize1, 'HorizontalAlignment', 'center');
            end
        end
    end

    plot((1) * (2 * T + 1) * 1.1 + [0, 0], -dh(1) + [0, yLineHight], 'Color', FrontColor, 'LineWidth', 1);
    plot((1) * (2 * T + 1) * 1.1 - [-T/3, +T/3], -dh(1) + [0, 0], 'Color', FrontColor, 'LineWidth', 1);
    plot((1) * (2 * T + 1) * 1.1 - [-T/3, +T/3], -dh(1) + [yLineHight, yLineHight], 'Color', FrontColor, 'LineWidth', 1);
    text((1) * (2 * T + 1) * 1.15 + T/3, -dh(1) + yLineHight / 2, [num2str(yLineHight) ''], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'Color', FrontColor, 'FontSize', Fontssize1);

    axis tight;
    hAx = get(h, 'CurrentAxes');
    set(hAx, 'Color', BackColor);
    set(hAx, 'XColor', BackColor);
    set(hAx, 'YColor', BackColor);
    set(hAx, 'ZColor', BackColor);
    set(hAx, 'XTickLabel', []);
    set(hAx, 'YTickLabel', []);
    title(['MU ' MU_ID], 'Color', FrontColor, 'FontSize', Fontssize1);
end