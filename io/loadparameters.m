[param, paramvalue, description] = textread(options.paramsfile, '%s %n %s', 'headerlines', 1, 'delimiter', '\t');
for li = 1 : length(param)
	str = [ param{li} ' = ' num2str(paramvalue(li)) ';' ];
	eval(str);
end

[copy_number, r_levels, description] = textread(options.levelsfile, '%n %n %s', 'headerlines', 1, 'delimiter', '\t');
for i = 1 : length(copy_number)
	LRR_levels(copy_number(i)+1) = r_levels(i);
	m_r0(copy_number(i)+1) = r_levels(i);
end

r_std = mad(r);

d 			= 2;
nStates     = 6;

z{1} = [ 1 ];
z{2} = [ 1 3 ];
z{3} = [ 1 2 3 ];
z{4} = [ 1 3 ];
z{5} = [ 1 2 2 3 ];
z{6} = [ 1 2 2 2 3 ];

cnloh = [ 0 0 0 1 0 0 0 ];

gnLabel{1} = 'NN';
gnLabel{2} = { 'A' 'B' 'NC' };
gnLabel{3} = { 'AA' 'AB' 'BB' 'NC' };
gnLabel{4} = { 'AA' 'BB' 'NC' };
gnLabel{5} = { 'AAA' 'AAB' 'ABB' 'BBB' 'NC' };
gnLabel{6} = { 'AAAA' 'AAAB' 'AABB' 'ABBB' 'BBBB' 'NC' };

m_b{1} = [ 0.5 ];
m_b{2} = [ 0 1 ];
m_b{3} = [ 0 1/2 1 ];
m_b{4} = [ 0 1 ];
m_b{5} = [ 0 1/3 2/3 1 ];
m_b{6} = [ 0 1/4 1/2 3/4 1 ];

w0{1} = 1;                          % homozygous deletion
w0{2} = [ 1/2 1/2 ];                % hemizygous deletion
w0{3} = [ 1/3 1/3 1/3 ];            % normal
w0{4} = [ 1/2 1/2 ];                % copy neutral loh
w0{5} = [ 1/3 1/6 1/6 1/3 ];        % duplication
w0{6} = [ 1/3 1/9 1/9 1/9 1/3 ];    % 4n balanced duplication

q0{1} 	= [ ones(1, nComp) ]./nComp;
q0{2} 	= [ ones(1, nComp) ]./nComp;
q0{3} 	= [ ones(1, nComp) ]./nComp;

d_r0{1} = zeros(1, nComp);
d_r0{2} = zeros(1, nComp);
d_r0{3} = zeros(1, nComp);

d_b0{1} = zeros(1, nComp);
d_b0{2} = zeros(1, nComp);
d_b0{3} = zeros(1, nComp);

S0_r_homdel 	= max(var( r(r < -1) ), r_std^2 );
S0_r_cnv 		= r_std^2;

S0{1}(:, :, 1) = [ 1e-6^2 0; 0 r_std^2 ];
S0{1}(:, :, 2) = [ 0.01^2 0; 0 r_std^2 ];
S0{1}(:, :, 3) = [ 0.05^2 0; 0 r_std^2 ];

S0{2}(:, :, 1) = [ 0.03^2 0; 0 r_std^2 ];
S0{2}(:, :, 2) = [ 0.04^2 0; 0 r_std^2 ];
S0{2}(:, :, 3) = [ 0.05^2 0; 0 r_std^2 ];

S0{3}(:, :, 1) = [ 1e-6^2 0; 0 r_std^2 ];
S0{3}(:, :, 2) = [ 0.01^2 0; 0 r_std^2 ];
S0{3}(:, :, 3) = [ 0.05^2 0; 0 r_std^2 ];

prior = ones(nStates+1, 1)/(nStates+1);

w 			= w0;
nu 			= 0.001;
q   		= q0;
d_r 		= d_r0;
d_b 		= d_b0;
S 			= S0;
S_r_homdel  = S0_r_homdel;
S_r_cnv     = S0_r_cnv;

options.nComp = nComp;
options.nStates = nStates;

initparams.nu = nu;
initparams.w = w;
initparams.q = q;
initparams.d_r = d_r;
initparams.d_b = d_b;
initparams.S = S;
initparams.S_r_homdel = S_r_homdel;
initparams.S_r_cnv = S_r_cnv;
initparams.L = options.L;
initparams.m_r = m_r0;

hyperparams.prior = prior;
hyperparams.w_alpha = w_alpha;
hyperparams.w0 = w0;
hyperparams.q_alpha = q_alpha;
hyperparams.q0 = q0;
hyperparams.nu_alpha = nu_alpha;
hyperparams.nu_beta = nu_beta;
hyperparams.tau = tau;
hyperparams.d_r0 = d_r0;
hyperparams.d_b0 = d_b0;
hyperparams.m_r0 = m_r0;
hyperparams.m_b = m_b;
hyperparams.S_alpha = S_alpha;
hyperparams.S_alpha_homdel = S_alpha_homdel;
hyperparams.S0 = S0;
hyperparams.S0_r_cnv = S0_r_cnv;
hyperparams.S0_r_homdel = S0_r_homdel;
hyperparams.v = v;
hyperparams.d = d;
hyperparams.z = z;
hyperparams.copy_number = copy_number;
hyperparams.r_max = r_max;
hyperparams.r_min = r_min;
hyperparams.gnLabel = gnLabel;
hyperparams.maxCopy = 5;
hyperparams.nComp = nComp;
hyperparams.nStates = nStates;
hyperparams.cnloh = cnloh;
hyperparams.longChromosome = longChromosome;