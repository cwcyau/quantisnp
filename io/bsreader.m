function bsreader(options)

headerlines = 10;

if ~isempty(options.genderfile)
	disp('QuantiSNP: Reading gender file');
	[sampleid, gender] = textread(options.genderfile, '%s\t%s', 'headerlines', 1);
end

% read beadStudio report file
disp(['QuantiSNP: Reading BeadStudio report header information']);
fid = fopen(options.bsfile, 'r');
columnOrder = zeros(1, 6);
columnNames = { 'Sample ID', 'SNP Name', 'Chr', 'Position', 'Log R Ratio', 'B Allele Freq' };
columnTypes = { ' %s' ' %s' ' %s' ' %f' ' %f' ' %f' };
for hi = 1 : headerlines
    
    tline = fgetl(fid);
    if strfind(tline, 'Num SNPs')
        disp(tline);
        [tmp, nSnps] = strread(tline, '%s\t%d', 'delimiter', '\t');
    end
    if strfind(tline, 'Num Samples')
        disp(tline);
        [tmp, nSamples] = strread(tline, '%s\t%d', 'delimiter', '\t');
    end
    if strfind(tline, 'Total SNPs')
        disp(tline);
        [tmp, totalSnps] = strread(tline, '%s\t%d', 'delimiter', '\t');
    end   
    if strfind(tline, 'SNP Name')
        disp(tline);
        str = strread(tline, '%s', 'delimiter', '\t');
        linedat = regexp(str, '\t', 'split');
        nColumns = length(linedat);
		format_str = [];
        for ci = 1 : nColumns
			matched = 0;
            for t = 1 : 6
                if strmatch(columnNames{t}, linedat{ci})
                    columnOrder(t) = ci;
					format_str = [ format_str columnTypes{t} ];
					matched = 1;
                end            
            end
			if matched == 0
				format_str = [ format_str ' %*s' ];
			end
        end
    end   
    
end
totalSnps = nSamples*nSnps;

%% summarise file data
disp([ 'QuantiSNP: Number of SNPs per sample: ' num2str(nSnps) ]);
disp([ 'QuantiSNP: Number of samples: ' num2str(nSamples) ]);
disp([ 'QuantiSNP: Total number of SNPs: ' num2str(totalSnps) ]);
disp('QuantiSNP: Variable - Column Number');
for ci = 1 : 6
    disp([ 'QuantiSNP: ' columnNames{ci} ' - ' num2str(columnOrder(ci)) ]);
end

if length(find(columnOrder > 0)) < 6
    disp('QuantiSNP: Error reading input file.');
    loc = find( columnOrder == 0 );
    for ci = 1 : length(loc)
        disp(['QuantiSNP: Check that the following column - ' columnNames{loc(ci)} ' - exists.']);
    end
    return;
end

%% for each sample

% start log file
try
	eid = fopen(options.logfile, 'wt');
catch
	disp(['QuantiSNP: Unable to write to: ' options.logfile]);
	return
end

fprintf(eid, 'SampleID\tProcessed\n');

	disp(['QuantiSNP: Reading BeadStudio data']);
	for si = 1 : nSamples

		% read data for each sample
		blockdat = textscan(fid, format_str, nSnps);  
		[tmp, I] = sort(columnOrder);
		sid = blockdat{I(1)};
		rs = blockdat{I(2)};
		chr = blockdat{I(3)};
		pos = blockdat{I(4)};
		r = blockdat{I(5)};
		b = blockdat{I(6)};

		% check if all the data comes from the same sample
		uid = unique(sid);
		if length(uid) > 1
			disp('QuantiSNP: Problem reading file. Input data format is incorrect. Check that the data file contains the same number of SNPs for every sample contained.');
			return;
		end
		disp(['QuantiSNP: Read data for sample: ' uid{1} ]);

		% find gender for this person
		options.sampleId = uid{1};
		options.gender = [];
		options.doGenderCalling = 0;
		if ~isempty(options.genderfile)
			loc = strmatch( uid{1}, sampleid, 'exact' );
			if isempty(loc)
				disp('QuantiSNP: Error! No gender for this sample was found. Using gender calling for this sample.');
				fprintf(eid, '%s\t%s\n', uid{1}, 'Yes - Gender Calling Applied');
				options.doGenderCalling = 1;
			else
				options.gender = gender{loc};
				disp([ 'QuantiSNP: Gender ' options.gender ' was specified for this sample ' uid{1} ]);
				switch lower(options.gender)
					case 'male'
						options.isMaleX = 1;
					case 'female'
						options.isMaleX = 0;
					otherwise
						disp('QuantiSNP: The gender specified must be male or female. Automatic gender calling will be used.');
						
				end					
				fprintf(eid, '%s\t%s\n', uid{1}, 'Yes');
				
			end
		else
			options.doGenderCalling = 1;
		end
		quantisnp2main(options, rs, chr, pos, r, b);

	end

	fclose(fid);

fclose(eid);

disp(['QuantiSNP: Log file written to: ' options.logfile ]);