function GN = genotypecalling(chrNo, vpath, genotype_alloc_prob, CNV, dataObj, hyperparams, options)

gnLabel = hyperparams.gnLabel;

GN = [];

rs = dataObj.rs;
pos = dataObj.pos;
b = dataObj.x(:, 1);
r = dataObj.x(:, 2);
snploc = dataObj.snploc;
cnvloc = dataObj.cnvloc;

nProbes = length(pos);

probeType = zeros(nProbes, 1, 'uint8');
probeType(cnvloc) = 1;

cnSeq = 2*ones(nProbes, 1, 'uint8');
if chrNo == options.chrX & options.isMaleX
	cnSeq(:) = 1;
end
bfSeq = NaN*zeros(nProbes, 1, 'uint8');

nCNV = length(CNV);
for cnvNo = 1 : nCNV
	index   = CNV{cnvNo}.index;
	bf      = CNV{cnvNo}.delta;
	cn      = CNV{cnvNo}.copy;
	cnSeq(index(1):index(2)) = cn;
	bfSeq(index(1):index(2)) = max(bf);
end

GN.rs 		= cell(nProbes, 1);
GN.cnSeq 	= zeros(nProbes, 1, 'uint8');
GN.bfSeq 	= zeros(nProbes, 1, 'uint8');
GN.ggn_txt 	= cell(nProbes, 1);
GN.ggn_prob = zeros(nProbes, 1, 'single');
GN.gn_txt 	= cell(nProbes, 1);
GN.gn_prob 	= zeros(nProbes, 1, 'single');

for t = 1 : nProbes
	
	if probeType(t) == 1
	
		ggn_txt 	= 'NA';
		gn_txt 		= 'NA';
		gn_prob 	= -1;
		ggn_prob 	= -1;
		
	else

		if vpath(t) ~= 1 & vpath(t) ~= hyperparams.nStates+1
			[ ggn_prob, ggn ] = max( genotype_alloc_prob{vpath(t)}(:, t) );
			ggn_txt = gnLabel{vpath(t)}{ggn};            
			[ gn_prob, gn ] = max( genotype_alloc_prob{3}(:, t) );	            
			gn_txt = gnLabel{3}{gn};
		else
			ggn_txt 	= 'NN';
			gn_txt 		= 'NN';
			gn_prob 	= -1;
			ggn_prob 	= -1;
		end
		
	end
	
	GN.rs{t}		= rs{t};
	GN.cnSeq(t) 	= cnSeq(t);
	GN.bfSeq(t) 	= bfSeq(t);
	GN.ggn_txt{t} 	= ggn_txt;
	GN.ggn_prob(t) 	= ggn_prob;
	GN.gn_txt{t} 	= gn_txt;
	GN.gn_prob(t) 	= gn_prob;
	
end
