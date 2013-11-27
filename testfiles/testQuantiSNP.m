 OUTDIR = '/data/suzaku/yau/luis/output/out/';
 INDIR = '/data/suzaku/yau/luis/output/quantisnp/';
 EMITERS = 10;
 LSETTING = 200000;
 PARAMSFILE = 'params.dat';
 CHRX = 23;
 CHRRANGE = '[1:22]';
 

 % ILLUMINA
P = genpath('/data/suzaku/yau/QuantiSNP2/release/03042011/source/');
addpath(P);
 
LEVELSFILE = 'levels-axiom.dat';
GCDIR = '/data/suzaku/yau/cnv/gc/b37/';

GENDER = 'female';


dirContents = dir(INDIR);

for fi = 162 : length(dirContents)

	fi

	INFILE = [ INDIR '/' dirContents(fi).name ];

	SAMPLEID = dirContents(fi).name(1:end-4);

	disp(SAMPLEID)

	quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');

end

break

 % ILLUMINA
P = genpath('/data/suzaku/yau/QuantiSNP2/release/03042011/source/');
addpath(P);
 
 LEVELSFILE = 'levels-hd.dat';
 GCDIR = '/data/suzaku/yau/cnv/gc/b36/';
 INFILE = [ '/data/suzaku/yau/downloads/sample.NA09748.res.nf' ];
 GENDER = 'female';
 
 SAMPLEID = '3p';
 quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
 
 

break


%  OUTDIR = '/data/suzaku/yau/QuantiSNP2/data/miroslava/';
%  EMITERS = 10;
%  LSETTING = 200000;
%  PARAMSFILE = 'params.dat';
%  CHRX = 23;
%  CHRRANGE = '[1:23]';
%  
%  % ILLUMINA
%  P = genpath('/data/suzaku/yau/QuantiSNP2/release/09042010/');
%  addpath(P);
%  
%  LEVELSFILE = 'levels-hd.dat';
%  GCDIR = '/data/suzaku/yau/cnv/gc/b36/';
%  INFILE = [ '/data/suzaku/yau/QuantiSNP2/data/miroslava/input 3p-.txt' ];
%  GENDER = 'female';
%  
%  SAMPLEID = '3p';
%  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  
%  
%  OUTDIR = '/data/suzaku/yau/QuantiSNP2/output/test/';
%  EMITERS = 10;
%  LSETTING = 200000;
%  PARAMSFILE = 'params.dat';
%  CHRX = 23;
%  CHRRANGE = '[1:23]';
%  
%  % ILLUMINA
%  P = genpath('/data/suzaku/yau/QuantiSNP2/release/09042010/');
%  addpath(P);
%  
%  %% quick test
%  CHRRANGE = '[6]';
%  
%  LEVELSFILE = 'levels.dat';
%  GCDIR = '/data/suzaku/yau/cnv/gc/b35/';
%  INFILE = [ '/data/suzaku/yau/cnv/WTCHG/HumanHap300/AA.txt' ];
%  GENDER = 'female';
%  
%  SAMPLEID = 'AA-with-gender-quick';
%  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  
%  
%  
%  
%  break
%  
%  %% test gender calling
%  CHRRANGE = '[1:23]';
%  
%  LEVELSFILE = 'levels.dat';
%  GCDIR = '/data/suzaku/yau/cnv/gc/b35/';
%  INFILE = [ '/data/suzaku/yau/cnv/WTCHG/HumanHap300/AA.txt' ];
%  GENDER = 'female';
%  
%  SAMPLEID = 'AA-with-gender';
%  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect', '--genotype');
%  
%  SAMPLEID = 'AA-no-gender';
%  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect', '--genotype');
%  %  
%  %  %% test gender calling
%  %  INFILE = [ '/data/suzaku/yau/cnv/WTCHG/HumanHap300/NM.txt' ];
%  %  GENDER = 'male';
%  %  
%  %  SAMPLEID = 'NM-with-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  SAMPLEID = 'NM-no-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  
%  %  %% test gender calling
%  %  INFILE = [ '/data/suzaku/yau/cnv/WTCHG/HumanHap300/OX1.txt' ];
%  %  GENDER = 'male';
%  %  
%  %  SAMPLEID = 'OX1-with-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  SAMPLEID = 'OX1-no-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  %% test gender calling
%  %  INFILE = [ '/data/suzaku/yau/cnv/WTCHG/HumanHap300/GM11840_1.txt' ];
%  %  GENDER = 'female';
%  %  
%  %  SAMPLEID = 'GM11840_1-with-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  SAMPLEID = 'GM11840_1-no-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  
%  %  %% test gender calling
%  %  INFILE = [ '/data/suzaku/yau/cnv/Scherer/1677734268.txt' ];
%  %  GENDER = 'female';
%  %  
%  %  SAMPLEID = '1677734268-with-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  SAMPLEID = '1677734268-no-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  
%  %  %% test Affy
%  %  LEVELSFILE = 'levels-affy.dat';
%  %  
%  %  GENDER = 'female';
%  %  INFILE = '/data/suzaku/yau/affy/output/hapmap/formatted/gw6.NA06985_GW6_C';
%  %  
%  %  SAMPLEID = 'NA06985-with-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--gender', GENDER, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--isaffy', '--chrX', CHRX, '--doXcorrect');
%  %  
%  %  SAMPLEID = 'NA06985-no-gender';
%  %  quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--isaffy', '--chrX', CHRX, '--doXcorrect');



OUTDIR = '/data/suzaku/yau/QuantiSNP2/output/test/';
EMITERS = 10;
LSETTING = 200000;
PARAMSFILE = 'params.dat';
CHRX = 23;
CHRRANGE = '[1:23]';

P = genpath('/data/suzaku/yau/QuantiSNP2/release/03042011/source/');
addpath(P);

%% test beadstudio reading
LEVELSFILE = 'levels-hd.dat';
GCDIR = '/data/suzaku/yau/cnv/gc/b36/';

%  GENDERFILE = '/data/suzaku/yau/QuantiSNP2/data/birgit/gender.txt';
%  BSFILE = '/data/suzaku/yau/QuantiSNP2/data/birgit/birgit_FinalReport.txt';
%  LOGFILE = '/data/suzaku/yau/QuantiSNP2/data/birgit/birgit.log';

%GENDERFILE = '/data/suzaku/yau/QuantiSNP2/data/pinto/gender.dat';
%BSFILE = '/data/suzaku/yau/QuantiSNP2/data/pinto/Test_FinalReport.txt';
%LOGFILE = '/data/suzaku/yau/QuantiSNP2/data/pinto/pinto.log';

BSFILE = '/data/suzaku/yau/downloads/TS_example.txt';
LOGFILE = '/data/suzaku/yau/downloads/TS_example.log';
OUTDIR = '/data/suzaku/yau/downloads/';

quantisnp2('--outdir', OUTDIR, '--logfile', LOGFILE, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--beadstudio-files', BSFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');

%  quantisnp2('--outdir', OUTDIR, '--logfile', LOGFILE, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--beadstudio-files', BSFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');
%  
%  BSFILE = '/data/suzaku/yau/QuantiSNP2/data/illuminaexample/HH550.txt';
%  LEVELSFILE = 'levels.dat';
%  GCDIR = '/data/suzaku/yau/cnv/gc/b35/';
%  LOGFILE = '/data/suzaku/yau/QuantiSNP2/output/test/HH550.log';
%  
%  quantisnp2('--outdir', OUTDIR, '--logfile', LOGFILE, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--beadstudio-files', BSFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');


