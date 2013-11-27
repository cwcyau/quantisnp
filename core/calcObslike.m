function [obslike_sum, genotype_alloc_prob, comp_alloc_prob, u] = calcObslike(dataObj, params, hyperparams)

if nargout < 3
	CALCCOMP = 0;
else
	CALCCOMP = 1;
end

x       	= dataObj.x;
snploc  	= dataObj.snploc;
cnvloc  	= dataObj.cnvloc;

nStates 	= hyperparams.nStates;
nComp 		= hyperparams.nComp;

m_b 		= hyperparams.m_b;
v 			= hyperparams.v;
z 			= hyperparams.z;
r_max 		= hyperparams.r_max;
r_min 		= hyperparams.r_min;

nu 			= params.nu;
w 			= params.w;
q 			= params.q;
m_r 		= params.m_r;
d_r 		= params.d_r;
d_b 		= params.d_b;
S 			= params.S;
S_r_homdel 	= params.S_r_homdel;
S_r_cnv 	= params.S_r_cnv;


n 			= size(x, 1);
nSnpProbes 	= length(snploc);
nCnvProbes 	= length(cnvloc);

obslike_sum 		= -1e100 + zeros(nStates+1, n );

genotype_alloc_prob = cell(nStates+1);
comp_alloc_prob 	= cell(nStates+1);

u 					= cell(nStates+1);
u{1}{1} 			= zeros(1, n);
u{nStates+1}{1} 	= zeros(1, n);
for i = 2 : nStates
	nGenotypes 					= length(w{i});
	genotype_alloc_prob{i} 		= zeros(nGenotypes+1, n);
	comp_alloc_prob{i} 			= cell(nGenotypes, 1);
	u{i} 						= cell(nGenotypes, 1);
	for j = 1 : nGenotypes
		comp_alloc_prob{i}{j} 	= zeros(nComp, n);
		u{i}{j} 				= zeros(nComp, n);
	end
end

[loglik, u{1}{1}] = logtpdf( x(:, 2), m_r(1), 1/S_r_homdel, v );
obslike_sum(1, :) = loglik;


for i = 2 : nStates
    
    nGenotypes = length(w{i});
	copy_number = hyperparams.copy_number(i);
   
    for j = 1 : nGenotypes

        mcType = z{i}(j);
        
        pr = zeros(nComp, n);
            
        for k = 1 : nComp
        
            iS  = inv(squeeze(S{mcType}(:, :, k)));

            if nSnpProbes > 0				
                m   = [ m_b{i}(j) + d_b{mcType}(k), m_r(copy_number+1) + d_r{mcType}(k)  ];
                [loglik, u{i}{j}(k, snploc)] = logtpdf(x(snploc, :), m, iS, v );
                pr(k, snploc) = log(1-nu) + log(w{i}(j)) + log(q{mcType}(k)) + loglik;                
            end
        
            if nCnvProbes > 0
                [loglik, u{i}{j}(k, cnvloc)] = logtpdf(x(cnvloc, 2), m_r(copy_number+1), 1/S_r_cnv, v );
                pr(k, cnvloc) = log(1-nu) + log(w{i}(j)) + log(q{mcType}(k)) + loglik;
            end 
            
        end
       
        pr_sum = logsumexp(pr, 1);
        
        obslike_sum(i, :) = logsumexp( [ obslike_sum(i, :); pr_sum ], 1);

        genotype_alloc_prob{i}(j, :) = pr_sum;

		if CALCCOMP
			comp_alloc_prob{i}{j} = pr;
			comp_alloc_prob{i}{j} = exp(comp_alloc_prob{i}{j} - repmat(pr_sum, [nComp 1]));
		end
     
    end
	
	% outlier state
    genotype_alloc_prob{i}(nGenotypes+1, :) = log(nu) - log(r_max-r_min);
    
    obslike_sum(i, :) = logsumexp( [ obslike_sum(i, :); genotype_alloc_prob{i}(nGenotypes+1, :) ], 1);
    
    gn_sum = logsumexp(genotype_alloc_prob{i}, 1);
    
    genotype_alloc_prob{i} = exp( genotype_alloc_prob{i} - repmat(gn_sum, [nGenotypes+1 1]) );
	
end

% amplification state
[loglik, u{nStates+1}{1}] = logtpdf( x(:, 2), m_r(hyperparams.maxCopy+1), 1/S_r_homdel, v  );
obslike_sum(nStates+1, :) = loglik;

