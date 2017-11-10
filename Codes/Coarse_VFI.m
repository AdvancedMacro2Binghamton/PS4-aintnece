clear all;
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

a_min=a_bar;
a_max=80;
num_a=50;
pas=(a_min+a_max)/num_a;
% Discretization of assets

a=a_min:pas:a_max;

% VFI
k_min=20;
k_max=40;

dis=1;

while abs(dis)>=0.01
    if abs(k_max-k_min)<1e-04
        break;
    end
    k_guess=(k_min+k_max)/2;
    r=alpha*(N_s/k_guess)^(1-alpha)+(1-delta);
    w=(1-alpha)*(k_guess/N_s)^alpha;
     % CURRENT RETURN (UTILITY) FUNCTION
   % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r*a', a);
    cons = bsxfun(@plus, cons, permute(z*w*l_bar, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(numel(z), numel(a));
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >1e-6
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:numel(z)
        E(:,:,i)=repmat(lzprob(i,:)*v_guess, [numel(a) 1]);
        end
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn, pol_indx]=max(v, [], 2);
        vfn=permute(vfn, [3 1 2]);
        v_tol = max(max(abs(vfn-v_guess)));
        v_guess=vfn;
    end
     % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
   % SET UP INITITAL DISTRIBUTION
    %Mu=ones(nz, num_a)/(nz*num_a);
    Mu = zeros(numel(z),numel(a));
  Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets
 % Mu = ones(nz, num_a); %alternative initial guess: same mass in all states
%Mu = Mu / sum(Mu(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1; 
while mu_tol > 1e-7
    [emp_ind, a_ind,mass] = find(Mu ); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (lzprob(emp_ind(ii), :)*Mu(emp_ind(ii), a_ind(ii)) )';
    end

    mu_tol = max(max(abs(MuNew-Mu)));
    
    Mu = MuNew ;
end
    aggsav=sum(Mu*a');
    dis=aggsav-k_guess;
    if dis>=0
       k_min=k_guess;
    else
       k_max=k_guess;
    end
    if abs(k_max-k_min)<10^(-5)
       break;
   end
end
% Euler equation error
c=bsxfun(@plus,bsxfun(@plus,-pol_fn,r*a),(z*w*l_bar)');
cf=c(:,pol_indx');
cf=reshape(cf,[numel(z) numel(a) numel(z)]);
for i=1:numel(z)
    c_prime(i,:)=lzprob(i,:)*cf(:,:,i);
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c_prime.^(-sigma)*r).*Mu));
