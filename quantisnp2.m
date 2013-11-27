function quantisnp2(varargin);

% --- display copyright notice --- %
str = sprintf('QuantiSNP v2.3\n----------------\n');
disp(str);
str = sprintf('\nThis software was developed and compiled using MATLAB R2008b \nand MATLAB Compiler 4.8 (R) (C) 1984-2010, The Mathworks, Inc.');
disp(str);
str = sprintf('\nBy using this software you agree to adhere to the terms and conditions\nof the licence supplied with this software.\n');
disp(str);
str = sprintf('\nCopyright (c) University of Oxford 2010. Website: http://www.well.ox.ac.uk/QuantiSNP/');
disp(str);
str = sprintf('\n-----------------------------------------------------------------------------\n');
disp(str);

% Common error messages
exitStr = ['QuantiSNP: Finished with an error.'];

nArgs = nargin;

chrRange		= 1:23;
isMaleX			= 0;
L				= 2e6;
EMiters 		= 10;
doGCcorrect 	= 0;
isAffy 			= 0;
doPlot			= 0;
doGenotyping	= 0;
doXcorrect 		= 0;
chrX			= 23;
doGenderCalling = 1;
doVerbose 		= 0;

paramsfile 		= 'params.dat';
levelsfile 		= 'levels.dat';
gcdir 			= [];
outdir 			= [];

bsfile 			= [];
logfile			= [];
genderfile 		= [];

sampleId		= [];
gender 			= [];
infile			= [];

for i = 1 : nargin

	%% required arguments
	if strmatch( '--outdir', lower(varargin{i}), 'exact')
		outdir = varargin{i+1};
	end
	
	%% optional arguments
	if strmatch( '--gcdir', lower(varargin{i}), 'exact')
		gcdir = varargin{i+1};
	end
	
	if strmatch( '--emiters', lower(varargin{i}), 'exact')
		EMiters = varargin{i+1};
	end

	if strmatch( '--lsetting', lower(varargin{i}), 'exact')
		L = varargin{i+1};
	end	
	
	if strmatch( '--isaffy', lower(varargin{i}), 'exact')
		isAffy = 1;
	end
	
	if strmatch('--chr', lower(varargin{i}), 'exact');
		chrRange = eval( [ '[' varargin{i+1} ']' ] );
	end
	
	if strmatch('--plot', lower(varargin{i}), 'exact');
		doPlot = 1;
	end

	if strmatch('--genotype', lower(varargin{i}), 'exact');
		doGenotyping = 1;
	end	
	
	if strmatch( '--config', lower(varargin{i}), 'exact')
		paramsfile = varargin{i+1};
	end

	if strmatch( '--levels', lower(varargin{i}), 'exact')
		levelsfile = varargin{i+1};
	end	
	
	if strmatch( '--doxcorrect', lower(varargin{i}), 'exact')
		doXcorrect = 1;
	end	
	
	if strmatch( '--chrx', lower(varargin{i}), 'exact')
		chrX = varargin{i+1};
	end	
	
	%% input specifications

	% single file processing
	if strmatch( '--input-files', lower(varargin{i}), 'exact')
		infile = varargin{i+1};
	end
	if strmatch( '--gender', lower(varargin{i}), 'exact')
		gender = varargin{i+1};
	end	
	if strmatch( '--sampleid', lower(varargin{i}), 'exact')
		sampleId = varargin{i+1};
	end

	% batch file processing
	if strmatch( '--beadstudio-files', lower(varargin{i}), 'exact')
		bsfile = varargin{i+1};
	end
	
	if strmatch( '--genderfile', lower(varargin{i}), 'exact')
		genderfile = varargin{i+1};
	end

	if strmatch( '--logfile', lower(varargin{i}), 'exact')
		logfile = varargin{i+1};
	end
	
	if strmatch( '--verbose', lower(varargin{i}), 'exact');
		doVerbose = 1;
	end

end

%% check output directory
if isempty(outdir)
	disp('QuantiSNP: No output directory specified.');
	return;
end

if ~exist(outdir, 'dir')
	disp(['QuantiSNP: The specified output directory ' outdir ' is invalid or does not exist.']);
	return;
end


%% check input specifications

% has any input file been specified
if isempty(infile) & isempty(bsfile)
	disp('QuantiSNP: No input file specified.');
	return;
end

if ~isempty(infile) & ~isempty(bsfile)
	disp('QuantiSNP: You may specify either single input file processing using --input-files or Beadstudio report processing using --beadstudio-files but not both.');
	return;
end

% which mode are we using?
doBSprocessing = 0;
if ~isempty(infile) & isempty(bsfile)
	doBSprocessing = 0;
	disp('QuantiSNP: Single-file mode input found.');
	disp(['QuantiSNP: Processing file: ' infile ]);
end

if isempty(infile) & ~isempty(bsfile)
	doBSprocessing = 1;
	disp('QuantiSNP: Beadstudio-file mode input found.');
	disp(['QuantiSNP: Processing file: ' bsfile ]);
end

% check the required arguments are supplied for each mode
if doBSprocessing == 1
	
	if isempty(logfile)
		disp(['QuantiSNP: Please specified a log file using the switch --logfile.']);
		return;
	end
	
	if isempty(genderfile)
		disp(['QuantiSNP: No gender file specified. Automatic gender calling will be used.']);
		doGenderCalling = 1;
	else
		doGenderCalling = 0;
	end
	
else

	if isempty(sampleId)
		disp(['QuantiSNP: Please specified a sample id using the switch --sampleid.']);
		return;
	end

	if isempty(gender)
		disp(['QuantiSNP: No gender specified. Automatic gender calling will be used.']);
		doGenderCalling = 1;
	else
		doGenderCalling = 0;
		switch lower(gender)
			case 'male'
				isMaleX = 1;
			case 'female'
				isMaleX = 0;
			otherwise
				disp('QuantiSNP: The gender specified must be male or female.');
				return;
		end		
	end	
	
end


%% check optional arguments

% check for local GC content files directory
if ~isempty(gcdir)
	disp('QuantiSNP: Local GC content directory specified. Local GC content correction will be used.');
	doGCcorrect = 1;
	if ~exist(gcdir, 'dir')
		disp([ 'QuantiSNP: The local GC content directory specified ' gcdir ' cannot be found or is invalid.' ]);
		return;
	end
else
	doGCcorrect = 0;
end

% check for parameters file
if ~isempty(paramsfile)
	if ~exist(paramsfile, 'file')
		disp(['QuantiSNP: The parameter file (' paramsfile ') was not found.']);
		return;
	end
end

% check for levels file
if ~isempty(levelsfile)
	if ~exist(levelsfile)
		disp(['QuantiSNP: The Log R Ratio levels file (' levelsfile ') was not found.']);
		return;
	end
end

% check that EMiters is numeric and greater than 0
str1 = ['QuantiSNP: The number of EM iterations must be a number greater than 0.'];
[Fail, EMiters] = checkArgs(EMiters, 'numeric', [1 Inf], str1, [], exitStr);
if Fail == 1
    return;
end

% check that Lsetting is numeric and greater than 0
str1 = ['QuantiSNP: The Lsetting must be a number greater than 0.'];
[Fail, L] = checkArgs(L, 'numeric', [1 Inf], str1, [], exitStr);
if Fail == 1
    return;
end

% check the ChrX number
str1 = ['QuantiSNP: The ChrX number must be greater than 0.'];
[Fail, chrX] = checkArgs(chrX, 'numeric', [1 Inf], str1, [], exitStr);
if Fail == 1
    return;
end


%% compile run-time options into one structure
options.doGenderCalling 	= doGenderCalling;
options.sampleId			= sampleId;
options.infile 				= infile;
options.isMaleX 			= isMaleX;
options.EMiters 			= EMiters;
options.doGCcorrect 		= doGCcorrect;
options.gcdir 				= gcdir;
options.chrRange 			= chrRange;
options.isAffy 				= isAffy;
options.L 					= L;
options.doPlot 				= doPlot;
options.doGenotyping 		= doGenotyping;
options.paramsfile 			= paramsfile;
options.levelsfile 			= levelsfile;
options.doXcorrect 			= doXcorrect;
options.chrX 				= chrX;
options.bsfile				= bsfile;
options.genderfile			= genderfile;
options.doBSread			= doBSprocessing;
options.logfile				= logfile;
options.outdir 				= outdir;
options.gender				= gender;
options.doVerbose			= doVerbose;

%% call main program
if doBSprocessing
	% do beadstudio file processing
	bsreader(options);
else
	% read file
	[rs, chr, pos, r, b] = textread(infile, '%s %s %n %n %n', 'headerlines', 1);
	
	% call quantisnp
	quantisnp2main(options, rs, chr, pos, r, b);
end

close all;
clear all;

disp('QuantiSNP: Exiting ...');