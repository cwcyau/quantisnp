function generateQC(params, options, hyperparams)

fid = fopen(options.outfile_qc, 'wt');
	fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'Sample ID', 'Chromosome', 'Outlier Rate', 'Std. Dev. LRR', 'Std. Dev. BAF', 'Gender');

for chrNo = options.chrRange

	nu = params{chrNo}.nu;

	nComp = length(params{chrNo}.d_r{2});

	n = 10000;

	rn = zeros(2, n);
	z = randsample([1:nComp], n, 'true', params{chrNo}.q{2});
	for i = 1 : n
		s = chol(squeeze(params{chrNo}.S{2}(:, :, z(i))));
		rn(:, i) = s*trnd(hyperparams.v, 2, 1);
		rn(1, i) = rn(1, i) + params{chrNo}.d_b{2}(z(i));
		rn(2, i) = rn(2, i) + params{chrNo}.d_r{2}(z(i));
	end	
	
	r_var = mad(rn(2, :));
	b_var = mad(rn(1, :));

	fprintf(fid, '%s\t%2.0f\t%1.4f\t%1.4f\t%1.4f\t%s\n', options.sampleId, chrNo, nu, r_var, b_var, options.gender);
	
end
fclose(fid);
