function PlotToFile2(dataObj, CNV, options)

chrX = options.chrX;

markSz = 2;
lineSz = 1;
fontSz = 8;
BFthresh = 30;

mplot = 4;
nplot = 4;

firstPlot = 1;
plotNo = 1;
plotRow = 1;

figure(1); clf;
set(gcf, 'Renderer', 'Painters');

orient landscape;
for k = 1 : mplot*nplot
    subplot(mplot, nplot, k);
	set(gca, 'Box', 'On', 'FontSize', fontSz);
end
hnd = suptitle(options.sampleId);
set(hnd, 'Interpreter', 'none');


fprintf('QuantiSNP. Plotting for chromosome: ');
for chrNo = options.chrRange

	if isempty(dataObj{chrNo})
		continue;
	end

	fprintf('%s ', num2str(chrNo));

    x = dataObj{chrNo}.x;
    
    pos = dataObj{chrNo}.pos/1e6;
    b = x(:, 1);
    r = x(:, 2);
    n = length(pos);
    
	if isempty(CNV{chrNo})
		continue;
	end    
    nCNV = length(CNV{chrNo});
    for i = 1 : nCNV
        
        ind = CNV{chrNo}{i}.index;
        cnvloc = CNV{chrNo}{i}.location;
        
        len = cnvloc(2) - cnvloc(1);
        nmarkers = ind(2) - ind(1) + 1;
        
        CN = CNV{chrNo}{i}.copy;
        BF = CNV{chrNo}{i}.delta(CNV{chrNo}{i}.type);
		if BF < 0 | ( CN == 2 & BF ~= 0 )
			continue;
		end

        markersOut = max(50, 2*nmarkers);
        
        cnv_range = ind(1):ind(2);
        r_cnv = r(cnv_range);
        b_cnv = b(cnv_range);
        pos_cnv = pos(cnv_range);
        
        loc = [ (ind(1)-markersOut):ind(1) ind(2):(ind(2)+markersOut) ];
        loc = loc(find(loc > 0 & loc <= n));
        r_out = r(loc);
        b_out = b(loc);
        pos_out = pos(loc);
        range_out = [ min(pos_out) max(pos_out) ];
        
        col = 'b';
        
        if CN < 2 & chrNo < chrX
            col = 'r';
        end
        if CN > 2 & chrNo < chrX
            col = 'g';
        end
        if CN < 1 & chrNo == chrX & options.isMaleX == 1
            col = 'r';
        end
        if CN > 1 & chrNo == chrX & options.isMaleX == 1
            col = 'g';
        end
        
        if CN < 2 & chrNo == chrX & options.isMaleX == 0
            col = 'r';
        end
        if CN > 2 & chrNo == chrX & options.isMaleX == 0
            col = 'g';
        end
        
        subplot(mplot, nplot, 2*(plotRow-1)*nplot + plotNo);
        hold on;
        plot(pos_out, r_out, 'ko', 'MarkerSize', markSz);
        plot(pos_cnv, r_cnv, 'o', 'Color', col, 'MarkerSize', markSz);
        xlim(range_out);
        ylim([-2 2]);
        if CN == 0
            ylim([-8 2]);
        end
        ylabel('LRR', 'FontSize', fontSz);
        set(gca, 'Box', 'On', 'FontSize', fontSz, 'XTick', linspace(range_out(1), range_out(2), 3));
        
        hnd1 = title(['Chr' num2str(chrNo) ':' num2str(cnvloc(1)) '-' num2str(cnvloc(2))]);
        set(hnd1, 'FontSize', fontSz);
        if BF > BFthresh
            set(hnd1, 'Color', 'r');
        end
        
        subplot(mplot, nplot, (2*(plotRow-1)+1)*nplot + plotNo);
        hold on;
        plot(pos_out, b_out, 'ko', 'MarkerSize', markSz);
        plot(pos_cnv, b_cnv, 'o', 'Color', col, 'MarkerSize', markSz);
        xlim(range_out);
        ylim([0 1]);
        xlabel('Position', 'FontSize', fontSz );
        ylabel('BAF', 'FontSize', fontSz);
        set(gca, 'Box', 'On', 'FontSize', fontSz, 'XTick', linspace(range_out(1), range_out(2), 3));
        
        hnd2 = title(['. CN: ' num2str(CN) '. Log BF: ' num2str(BF) ]);
        set(hnd2, 'FontSize', fontSz);
        if BF > BFthresh
            set(hnd2, 'Color', 'r');
        end

        plotNo = plotNo + 1;
        
        if plotNo > nplot
            plotNo = 1;
            plotRow = plotRow + 1;
        end
        
        if ( plotRow > 1/2*mplot )
            if firstPlot == 1
                print(gcf, options.outfile_plot, '-r300', '-dpsc2');
                firstPlot = 0;
				figure(1); clf;
            else
                print(gcf, options.outfile_plot, '-r300', '-dpsc2', '-append');
				figure(1); clf;
            end
            
            if ~( ( chrNo == options.chrRange(end) ) & ( i == nCNV ) )
                clf;
                for k = 1 : mplot*nplot
                    subplot(mplot, nplot, k);
                    set(gca, 'Box', 'On', 'FontSize', fontSz);
                end
                plotRow = 1;
            end
            
        end
        
    end
    
end
fprintf('\n');

if firstPlot == 1
    print(gcf, options.outfile_plot, '-r300', '-dpsc2');
    firstPlot = 0;
else
    if plotNo > 1
        print(gcf, options.outfile_plot, '-r300', '-dpsc2', '-append');
    end
end

close all;

gzip(options.outfile_plot);
delete(options.outfile_plot);
