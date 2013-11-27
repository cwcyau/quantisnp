function params = updateParams(expectationObj, dataObj, params, hyperparams, options)

chrRange        = options.chrRange;

nStates         = hyperparams.nStates;
nComp           = hyperparams.nComp;

z               = hyperparams.z;
d               = hyperparams.d;

nu_alpha        = hyperparams.nu_alpha;
nu_beta         = hyperparams.nu_beta;

w_alpha         = hyperparams.w_alpha;
w0              = hyperparams.w0;

q_alpha         = hyperparams.q_alpha;
q0              = hyperparams.q0;

tau             = hyperparams.tau;
m_r0             = hyperparams.m_r0;
m_b             = hyperparams.m_b;
d_r0            = hyperparams.d_r0;
d_b0            = hyperparams.d_b0;

S0              = hyperparams.S0;
S0_r_cnv        = hyperparams.S0_r_cnv;
S0_r_homdel     = hyperparams.S0_r_homdel;
S_alpha         = hyperparams.S_alpha;
S_alpha_homdel  = hyperparams.S_alpha_homdel;

m_r = params.m_r;

nu_num = 0;
nu_den = 0;

for i = 2 : nStates

	nGenotypes = length(hyperparams.z{i});

    w_num = zeros(1, nGenotypes);
    w_den = zeros(1, nGenotypes);
	
    for chrNo = chrRange

		if isempty(dataObj{chrNo}) 
			continue;
		end
        
        pos = dataObj{chrNo}.pos;
        x = dataObj{chrNo}.x;
		snploc = dataObj{chrNo}.snploc;
        
        gamma = expectationObj{chrNo}.gamma;
        genotype_alloc_prob = expectationObj{chrNo}.genotype_alloc_prob;
        
        nProbes = length(pos);
        
        if nProbes < 2
            continue;
        end
        
        for j = 1 : nGenotypes+1
            
            pr = gamma(i, :).*genotype_alloc_prob{i}(j, :);
            gn_sum = sum(pr);
            
            if j == nGenotypes+1
                nu_num = nu_num + gn_sum;
            end
            nu_den = nu_den + gn_sum;

        end
        
        for j = 1 : nGenotypes
            
            pr = gamma(i, :).*genotype_alloc_prob{i}(j, :);
            gn_sum = sum(pr);

			w_num(j) = w_num(j) + gn_sum;
            
        end
        
    end
    
    w{i} = (w_num + w_alpha*w0{i} - 1)./( sum(w_num) + sum(w_alpha*w0{i}) - nGenotypes );
	
end

nu = (nu_num + nu_alpha)/(nu_den + nu_beta);


for mcType = 1 : 3
    q_num{mcType}   = zeros(1, nComp);
    d_r_num{mcType} = zeros(1, nComp);
    d_r_den{mcType} = zeros(1, nComp);
    d_b_num{mcType} = zeros(1, nComp);
    d_b_den{mcType} = zeros(1, nComp);
    S_num{mcType}   = zeros(d, d, nComp);
    S_den{mcType}   = zeros(d, d, nComp);
	d_num{mcType} 	= zeros(2, nComp);
	d_den{mcType} 	= zeros(nComp);
end

for chrNo = chrRange

	if isempty(dataObj{chrNo})
		continue;
	end
    
    x                   = dataObj{chrNo}.x;
    snploc              = dataObj{chrNo}.snploc;
    
    gamma               = expectationObj{chrNo}.gamma;
    genotype_alloc_prob = expectationObj{chrNo}.genotype_alloc_prob;
    comp_alloc_prob     = expectationObj{chrNo}.comp_alloc_prob;
    u                   = expectationObj{chrNo}.u;
    
    for i = 2 : nStates
    
		cn = hyperparams.copy_number(i);
        nGenotypes = length(z{i});
        
        for j = 1 : nGenotypes
            
            mcType = z{i}(j);
            
            for k = 1 : nComp
                
                pr = gamma(i, :).*genotype_alloc_prob{i}(j, :).*comp_alloc_prob{i}{j}(k, :);
                
                comp_sum = sum(pr);
                
                q_num{mcType}(k)   = q_num{mcType}(k) + comp_sum;
                
                pr = gamma(i, :).*genotype_alloc_prob{i}(j, :).*comp_alloc_prob{i}{j}(k, :).*u{i}{j}(k, :);

				if k > 1
				
					y = x(snploc, :)';
					y(1, :) = y(1, :) - m_b{i}(j);
					y(2, :) = y(2, :) - m_r(cn+1);

					d_num{mcType}(:, k) = d_num{mcType}(:, k) + sum( y.*repmat(pr(snploc), [2 1]) , 2);
					d_den{mcType}(k) = d_den{mcType}(k) + sum(pr(snploc));
				
				end

            end
            
        end
        
    end
        
end

for mcType = 1 : 3

    q{mcType} = (q_num{mcType} + q_alpha*q0{mcType} - 1)./(sum(q_num{mcType}) + sum(q_alpha*q0{mcType}) - nComp);
	for k = 2 : nComp
		dd_0 = [ d_b0{mcType}(k); d_r0{mcType}(k) ];
		dd = 1/(d_den{mcType}(k) + tau)*( d_num{mcType}(:, k) + tau.*dd_0 );
		d_b{mcType}(k) = dd(1);
		d_r{mcType}(k) = dd(2);
	end
	
end


%% update variances
S_r_homdel_num = 0;
S_r_homdel_den = 0;

S_r_cnv_num = 0;
S_r_cnv_den = 0;

for chrNo = chrRange

	if isempty(dataObj{chrNo})
		continue;
	end
   
    x                   = dataObj{chrNo}.x;
    snploc              = dataObj{chrNo}.snploc;
    cnvloc              = dataObj{chrNo}.cnvloc;
    
    gamma               = expectationObj{chrNo}.gamma;
    genotype_alloc_prob = expectationObj{chrNo}.genotype_alloc_prob;
    comp_alloc_prob     = expectationObj{chrNo}.comp_alloc_prob;
    u                   = expectationObj{chrNo}.u;
    
    for i = [ 1 nStates+1 ]
        
		cn = hyperparams.copy_number(i);
		
        pr = gamma(i, :).*u{i}{1};
        
        rd = x(:, 2)' - m_r(cn+1);
        
        S_r_homdel_num = S_r_homdel_num + sum(pr.*rd.^2);
        S_r_homdel_den = S_r_homdel_den + sum(pr);
        
    end
    
    for i = 2 : nStates
    
		cn = hyperparams.copy_number(i);
        nGenotypes = length(z{i});
        
        for j = 1 : nGenotypes
            
            mcType = z{i}(j);
            
            for k = 1 : nComp
                
                pr = gamma(i, :).*genotype_alloc_prob{i}(j, :).*comp_alloc_prob{i}{j}(k, :).*u{i}{j}(k, :);
                
                % for SNP probes
                rd = x(:, 2)' - m_r(cn+1) - d_r{mcType}(k);
                bd = x(:, 1)' - m_b{i}(j) - d_b{mcType}(k);
                
                S_num{mcType}(1, 1, k) = S_num{mcType}(1, 1, k) + sum(pr(snploc).*bd(snploc).^2);
                S_num{mcType}(1, 2, k) = S_num{mcType}(1, 2, k) + sum(pr(snploc).*bd(snploc).*rd(snploc));
                S_num{mcType}(2, 1, k) = S_num{mcType}(2, 1, k) + sum(pr(snploc).*bd(snploc).*rd(snploc));
                S_num{mcType}(2, 2, k) = S_num{mcType}(2, 2, k) + sum(pr(snploc).*rd(snploc).^2);
                
                S_den{mcType}(1, 1, k) = S_den{mcType}(1, 1, k) + sum(pr(snploc));
                S_den{mcType}(1, 2, k) = S_den{mcType}(1, 2, k) + sum(pr(snploc));
                S_den{mcType}(2, 1, k) = S_den{mcType}(2, 1, k) + sum(pr(snploc));
                S_den{mcType}(2, 2, k) = S_den{mcType}(2, 2, k) + sum(pr(snploc));
                
                % for CNV probes
                if length(cnvloc) > 0
                    S_r_cnv_num = S_r_cnv_num + sum(pr(cnvloc).*rd(cnvloc).^2);
                    S_r_cnv_den = S_r_cnv_den + sum(pr(cnvloc));
                end
                
            end
            
        end
        
    end
      
end

for mcType = 1 : 3

    for k = 1 : nComp
	
        S{mcType}(1, 1, k) = ( S_num{mcType}(1, 1, k) + tau*(d_b{mcType}(k)-d_b0{mcType}(k))^2 + S_alpha*S0{mcType}(1, 1, k) )./( S_den{mcType}(1, 1, k) + S_alpha - d );
		
        S{mcType}(1, 2, k) = ( S_num{mcType}(1, 2, k) + tau*(d_r{mcType}(k)-d_r0{mcType}(k))*(d_b{mcType}(k)-d_b0{mcType}(k)) + S_alpha*S0{mcType}(1, 2, k) )./( S_den{mcType}(1, 2, k) + S_alpha - d );
		
        S{mcType}(2, 1, k) = ( S_num{mcType}(2, 1, k) + tau*(d_r{mcType}(k)-d_r0{mcType}(k))*(d_b{mcType}(k)-d_b0{mcType}(k)) + S_alpha*S0{mcType}(2, 1, k) )./( S_den{mcType}(2, 1, k) + S_alpha - d );
		
        S{mcType}(2, 2, k) = ( S_num{mcType}(2, 2, k) + tau*(d_r{mcType}(k)-d_r0{mcType}(k))^2 + S_alpha*S0{mcType}(2, 2, k) )./( S_den{mcType}(2, 2, k) + S_alpha - d );
		
    end
end
S_r_homdel  = ( S_r_homdel_num + S_alpha_homdel*S0_r_homdel )/( S_r_homdel_den + S_alpha_homdel - 1 );
S_r_cnv     = ( S_r_cnv_num + S_alpha*S0_r_cnv )/( S_r_cnv_den + S_alpha - 1 );


%% update CN levels
m_r_num = zeros(1, max(hyperparams.copy_number)+1);
m_r_den = zeros(1, max(hyperparams.copy_number)+1);
for chrNo = chrRange

	if isempty(dataObj{chrNo})
		continue;
	end

	x                   = dataObj{chrNo}.x;
	snploc              = dataObj{chrNo}.snploc;

	gamma               = expectationObj{chrNo}.gamma;
	genotype_alloc_prob = expectationObj{chrNo}.genotype_alloc_prob;
	comp_alloc_prob     = expectationObj{chrNo}.comp_alloc_prob;
	u                   = expectationObj{chrNo}.u;

	for i = [ 1 : nStates+1 ]

		cn = hyperparams.copy_number(i);
		
		if cn == 2
			continue;
		end
		
		if cn == 0 | cn == hyperparams.maxCopy
		
			pr = gamma(i, :).*u{i}{1};

			m_r_num(cn+1) = m_r_num(cn+1) + sum(pr.*x(:, 2)');
			m_r_den(cn+1) = m_r_den(cn+1) + sum(pr);
		
		else
			
			nGenotypes = length(z{i});

			for j = 1 : nGenotypes

				mcType = z{i}(j);

				for k = 1 : nComp
				
					pr = gamma(i, :).*genotype_alloc_prob{i}(j, :).*comp_alloc_prob{i}{j}(k, :).*u{i}{j}(k, :);

					dd = 2*S{mcType}(1, 1, k)*(x(:, 2)'-d_r{mcType}(k)) - (S{mcType}(2, 1, k)+S{mcType}(1, 2, k))*(x(:, 1)'-d_b{mcType}(k));
					m_r_num(cn+1) = m_r_num(cn+1) + sum(pr.*dd);
					m_r_den(cn+1) = m_r_den(cn+1) + sum(pr*2*S{mcType}(1, 1, k));

				end

			end
		
		end

	end

end

m_r = ( m_r_num + tau*m_r0 )./( m_r_den + tau );


%% store
params.w            = w;
params.nu           = nu;
params.q            = q;
params.d_r          = d_r;
params.d_b          = d_b;
params.S            = S;
params.S_r_homdel   = S_r_homdel;
params.S_r_cnv      = S_r_cnv;
params.m_r      	= m_r;

