% Parameter setups
alpha=1/3;
beta=0.99;
sigma=2;
delta=0.025;
a_bar=0;
l_bar=1;
rho=0.5;
sigma_e=0.2;
sigma_z=sigma_e/sqrt(1-rho^2);

% Discretization using Tauchen's method
[lz,lzprob]=TAUCHEN(5,rho,sigma_e,3);
z=exp(lz');
P=lzprob^1000;
N_s=sum(P(1,:).*z);

% Discretization of assets
a_min=a_bar;
a_max=7;
num_a=12;

a=linspace(a_min,a_max,num_a);

% VFI
k_min=a_bar;
k_max=a_max;

dis=1;

    k_guess=(k_min+k_max)/2;
    r=alpha*(N_s/k_guess)^(1-alpha)+(1-delta);
    w=(1-alpha)*(k_guess/N_s)^alpha;
     % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus,r*a,a');
    cons = bsxfun(@plus, cons, permute(z*w*l_bar, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(numel(z), num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >1e-06
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        value_mat=ret+beta*permute(lzprob*v_guess,[3 2 1]);
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn,pol_indx]=max(value_mat,[],2);
        vfn=permute(vfn,[3 1 2]);
        v_tol=max(max(abs(vfn-v_guess)));
        v_guess=vfn;
  
    end
    
    % KEEP DECSISION RULE
    pol_indx = permute(pol_indx,[3 1 2]);
    pol_fn=a(pol_indx);