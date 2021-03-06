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
num_a=12;
a=linspace(a_min,a_max,num_a);

% VFI
k_min=20;
k_max=40;

dis=1;

dis=1;
aggsav=0;
while abs(dis)>=0.01
    if abs(k_max-k_min)<1e-04
        break;
    end
    k_guess=(k_min+k_max)/2;
    r=alpha*(N_s/k_guess)^(1-alpha)+(1-delta);
    w=(1-alpha)*(k_guess/N_s)^alpha;
   % CURRENT RETURN (UTILITY) FUNCTION
    cons_c = bsxfun(@minus, r*a', a);
    cons_c = bsxfun(@plus, cons_c, permute(z*w*l_bar, [1 3 2]));
    ret_c = (cons_c .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret_c(cons_c<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess_c = zeros(numel(z), numel(a));
    % VALUE FUNCTION ITERATION
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:numel(z)
        E_c(:,:,i)=repmat(lzprob(i,:)*v_guess_c, [numel(a) 1]);
        end
        value_mat_c=ret_c+beta*E_c;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn_c, pol_indx]=max(value_mat_c, [], 2);
        vfn_c=permute(vfn_c, [3 1 2]);
        a_interp=linspace(a_min,a_max,500);
        v_interp=interp1(a,vfn_c',a_interp);
        cons_interp = bsxfun(@minus, r*a_interp', a_interp);
        cons_interp = bsxfun(@plus, cons_interp, permute(z*w*l_bar, [1 3 2]));
        ret_interp = (cons_interp .^ (1-sigma)) ./ (1 - sigma);
        ret_interp(cons_interp<0)=-Inf;
        v_guess=v_interp';
        v_tol=1;
        while v_tol>1e-6
        for i=1:numel(z)
        E(:,:,i)=repmat(lzprob(i,:)*v_guess, [numel(a_interp) 1]);
        end
        value_mat=ret_interp+beta*E;
        [vfn,pol_indx]=max(value_mat,[],2);
        vfn=permute(vfn,[3 1 2]);
        v_tol = max(max(abs(vfn-v_guess)));
        v_guess=vfn;
        end
   pol_indx=permute(pol_indx, [3 1 2]);
   pol_fn=a_interp(pol_indx);
   % SET UP INITITAL DISTRIBUTION
   Mu = zeros(numel(z),numel(a_interp));
  Mu(1, 4) = 1;
  % ITERATE OVER DISTRIBUTIONS
    mu_tol=1;
    while mu_tol >1e-7
          [emp_ind, a_ind] = find(Mu); % find non-zero indices
          MuNew = zeros(size(Mu));
          for ii = 1:length(emp_ind)
              apr_ind = pol_indx(emp_ind(ii), a_ind(ii));
              MuNew(:, apr_ind) = MuNew(:, apr_ind) + (lzprob(emp_ind(ii),:)*Mu(emp_ind(ii),a_ind(ii)))';
          end
          mu_tol=max(max(abs(MuNew-Mu)));
          Mu=MuNew;
    end
    aggsav=sum(Mu*a_interp');
    dis=aggsav-k_guess;
    if dis>=0
       k_min=k_guess;
    else
       k_max=k_guess;
    end
    if abs(k_max-k_min)<1e-5
       break;
   end
end
% Euler equation error
c=bsxfun(@plus,bsxfun(@plus,-pol_fn,r*a_interp),(z*w*l_bar)');
cf=c(:,pol_indx');
cf=reshape(cf,[numel(z) numel(a_interp) numel(z)]);
for i=1:numel(z)
    c_prime(i,:)=lzprob(i,:)*cf(:,:,i);
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c_prime.^(-sigma)*r).*Mu));