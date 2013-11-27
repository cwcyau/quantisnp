function genotype(GN, dataObj, options, hyperparams)

gnLabel = hyperparams.gnLabel;

fid = fopen(options.outfile_genotype, 'wt');
	fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
				'ProbeName', 'Chr', 'Position', 'Log R Ratio', 'B Allele Frequency', ...
				'Copy Number', 'Max. BF', 'Generalised Genotype', 'Generalised Genotype Prob', ...
				'Diploid Genotype', 'Diploid Genotype Prob' );

for chrNo = options.chrRange
    
	if isempty(dataObj{chrNo})
		continue;
	end
	
    rs = dataObj{chrNo}.rs;
    pos = dataObj{chrNo}.pos;
    b = dataObj{chrNo}.x(:, 1);
    r = dataObj{chrNo}.x(:, 2);
    
    nProbes = length(pos);

    for t = 1 : nProbes
		
		try
			fprintf(fid,'%s\t%2.0f\t%15.0f\t%1.3f\t%1.3f\t%1.0f\t%10.3f\t%s\t%1.3f\t%s\t%1.3f\n', ...
				rs{t}, ...
				chrNo, ...
				pos(t), ...
				r(t), ...
				b(t), ...
				GN{chrNo}.cnSeq(t), ...
				GN{chrNo}.bfSeq(t), ...
				GN{chrNo}.ggn_txt{t}, ...
				GN{chrNo}.ggn_prob(t), ...
				GN{chrNo}.gn_txt{t}, ...
				GN{chrNo}.gn_prob(t) );
		catch
			chrNo
			t
			GN{chrNo}
		end
        
    end
    
end

fclose(fid);

%  gzip(options.outfile_genotype);
%  delete(options.outfile_genotype);
