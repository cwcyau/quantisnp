function [ dataObj, options ] = constructData(rs, chr, pos, r, b, hyperparams, options)

if options.doGCcorrect
	disp('QuantiSNP. Using local GC correction.');
end

if options.doXcorrect
	disp('QuantiSNP. Using ChrX correction.');
end

disp(['QuantiSNP. Chr' num2str(options.chrX) ' is the X chromosome' ]);

fprintf('QuantiSNP. Reading data for chromosome: ');
for chrNo = options.chrRange

	fprintf('%s ', num2str(chrNo));
	dataObj{chrNo} = [];	

    if chrNo ~= options.chrX
        chrloc = strmatch(num2str(chrNo), chr, 'exact');
        gcfile = [ options.gcdir '/' num2str(chrNo) '_1k.txt' ];
    else
        chrloc = strmatch('X', chr, 'exact');
        gcfile = [ options.gcdir '/' 'X_1k.txt' ];
    end

    if length(chrloc) < 10
        continue;
    end

	n_chr = length(chrloc);

	rs_chr = cell(n_chr, 1);
	pos_chr = zeros(n_chr, 1);
	x_chr = zeros(n_chr, 2);

    rs_chr  		= rs(chrloc);
    pos_chr 		= pos(chrloc);
    x_chr(:, 2) 	= r(chrloc);
    x_chr(:, 1) 	= b(chrloc);

    [pos_chr, I] 	= sort(pos_chr);
    rs_chr 			= rs_chr(I);
    x_chr(:, 2) 	= x_chr(I, 2);
    x_chr(:, 1) 	= x_chr(I, 1);

	snpprobes = ones(1, n_chr);
    if options.isAffy == 1
        cnvloc = strmatch('CN', rs_chr);
    else
        cnvloc = strmatch('cnv', rs_chr);
    end
    snpprobes(cnvloc) = 0;
    snploc = find( snpprobes == 1 );
	
    if length(cnvloc) > 1
        cnvProbesFound = 1;
    else
        cnvProbesFound = 0;
    end

    % do gc correction
    if options.doGCcorrect == 1
	
    	try
	        [gcpos, gcend, gc] = textread(gcfile, '%n %n %n');
			gc = gc/100;
			gc = gc - mean(gc);
        catch
        	disp(['QuantiSNP. Error! Unable to read GC content file. Stopping.']);
        	return;
        end
		
        gc_s = nanmoving_average(gc, 100);
        gc_chr = interp1(gcpos, gc_s, pos_chr, 'nearest', 'extrap');
        betas = robustfit(gc_chr, x_chr(:, 2));

        x_chr(:, 2) = x_chr(:, 2) - betas(2)*gc_chr;

    end

	% do gender calling
	if chrNo == options.chrX & options.doGenderCalling == 1
	
		nHetX = length( find( x_chr(:, 1) > 0.25 & x_chr(:, 1) < 0.75 ) );
		nX = length( x_chr(:, 1) );

		if nHetX/nX < 0.05
			options.gender = 'male';
			options.isMaleX = 1;
		else
			options.gender = 'female';
			options.isMaleX = 0;
		end
	end

    % do chrX correction
    if chrNo == options.chrX & options.doXcorrect == 1
        x_chr(:, 2) = x_chr(:, 2) - median(x_chr(:, 2));
        if options.isMaleX == 1
            x_chr(:, 2) = x_chr(:, 2) + hyperparams.m_r0(2);
        end
    end

    %% generate sub-samples by sampling from uniformly spaced quantiles
	[r_sorted, I] = sort(x_chr(:, 2));

	n = length(pos_chr);
	loc = [ 1 : n ];

	loc = loc(I);
	if options.doSubSample
		loc = loc(1:options.subSampleLevel:end);
	end
	loc = sort(loc);

	dataObj{chrNo}.rs    = rs_chr(loc);
	dataObj{chrNo}.pos   = pos_chr(loc);
	dataObj{chrNo}.x     = x_chr(loc, :);

	snpprobes = ones(1, n);
    if options.isAffy == 1
        dataObj{chrNo}.cnvloc = strmatch('CN', rs_chr(loc));        
    else
        dataObj{chrNo}.cnvloc = strmatch('cnv', rs_chr(loc));
    end
    snpprobes(dataObj{chrNo}.cnvloc) = 0;
    dataObj{chrNo}.snploc = find( snpprobes == 1 );


end
fprintf('\n');

if options.doGenderCalling == 1
	if options.isMaleX == 1
		disp('QuantiSNP: Found male sample');
	else
		disp('QuantiSNP: Found female sample');
	end
end
