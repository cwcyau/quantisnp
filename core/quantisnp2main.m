function quantisnp2main(options, rs, chr, pos, r, b)
%
% quantisnp2main(options, rs, chr, pos, r, b)
%

options.outfile = [ options.outdir '/' options.sampleId '.cnv'];
options.outfile_loh = [ options.outdir '/' options.sampleId '.loh'];
options.outfile_genotype = [ options.outdir '/' options.sampleId '.gn'];
options.outfile_plot = [ options.outdir '/' options.sampleId '.ps' ];
options.outfile_qc = [ options.outdir '/' options.sampleId '.qc' ];

startTime = cputime; 

% filter out bad data
loc = find( ~isnan(r) & ~isnan(b) & ~isinf(r) & ~isinf(b) );
rs = rs(loc);
chr = chr(loc);
pos = pos(loc);
r = r(loc);
b = b(loc);

% do median correction
loc = find( ~isnan(r) & ~isinf(b) );
r 	= r - median(r(loc));

% load hyperparameters
r_max = max(max(r), 1);
r_min = min(min(r), -4);

loadparameters;
paramsX = initparams;

%% process full dataset
chrRange = options.chrRange;
options.doSubSample = 0;
[ dataObj, options ] = constructData(rs, chr, pos, r, b, hyperparams, options);
clear rs chr pos r b;

chrRange

if length(chrRange) > 0

	params = cell(1, 100);
	for chrNo = chrRange
		fprintf('QuantiSNP. Using EM for parameter estimation. Chromosome: %g.\n', chrNo);
		options.chrRange = chrNo;
		params{chrNo} = learningEM(dataObj, initparams, hyperparams, options);
	end

end

%% do CNV calling
options.chrRange = chrRange;
[ CNV, GN ] = cnvcalling(dataObj, params, hyperparams, options);

%% generate QC file
fprintf('%s %s\n', 'QuantiSNP. Writing QC file: ', options.outfile_qc);
generateQC(params, options, hyperparams);

%% write to file
fprintf('%s %s\n', 'QuantiSNP. Writing output to file: ', options.outfile);
writeCNVtoFile(options, options.sampleId, chrRange, CNV);

%% write genotypes to file
if options.doGenotyping
	fprintf('%s %s\n', 'QuantiSNP. Writing genotypes to file: ', options.outfile_genotype);
	genotype(GN, dataObj, options, hyperparams);
end

%% plot to file
if options.doPlot
	fprintf('%s %s\n', 'QuantiSNP. Plotting to file: ', options.outfile_plot);
	PlotToFile2(dataObj, CNV, options);
end

endTime = cputime;

timeTaken = (endTime - startTime)/60;

str = [ 'QuantiSNP. Done in ' num2str(timeTaken, '%1.2f') ' mins.' ];
disp(str);
