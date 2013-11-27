 % ILLUMINA
P = genpath('/data/suzaku/yau/QuantiSNP2/release/03042011/source/');
addpath(P);

OUTDIR = '/data/suzaku/yau/temp/';
INDIR = '/data/suzaku/yau/luis/output/quantisnp/';
EMITERS = 1;
LSETTING = 200000;
PARAMSFILE = 'params.dat';
CHRX = 23;
CHRRANGE = '[1:26]';
 
 
LEVELSFILE = 'levels-axiom.dat';
GCDIR = '/data/suzaku/yau/cnv/gc/b37/';

GENDER = 'female';


dirContents = dir(INDIR);

for fi = 1 : length(dirContents)

	INFILE = [ INDIR '/' dirContents(fi).name ];

	SAMPLEID = dirContents(fi).name(1:end-4);

	disp(SAMPLEID)

	quantisnp2('--outdir', OUTDIR, '--sampleid', SAMPLEID, '--emiters', EMITERS, '--lsetting', LSETTING, '--gcdir', GCDIR, '--input-files', INFILE, '--chr', CHRRANGE, '--plot', '--levels', LEVELSFILE, '--params', PARAMSFILE, '--chrX', CHRX, '--doXcorrect');

end

