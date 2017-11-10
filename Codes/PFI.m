clear all;
l_bar=1;
sigma_e=0.2;
rho=0.5;
m=3;
nz=5;
[lz,lzprob] =TAUCHEN(nz,rho,sigma_e,m);
z=exp(lz');
P=lzprob^1000;
N_s=sum(P(1,:).*z);
a_min=a_bar;
a_max=80;
num_a=500;
pas=(a_min+a_max)/num_a;
% Discretization of assets

a=a_min:pas:a_max;

% VFI
k_min=20;
k_max=40;

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
        V=reshape(vfn, [numel(z)*numel(a) 1]);
        pol_indx1=permute(pol_indx, [3 1 2]);
        pol_fn1 = a(pol_indx1);
        cons1=r*repmat(a,[numel(z) 1])-pol_fn1;
        cons1 = cons1+repmat((z*w*l_bar)', [1 numel(a)]);
        ret1 = (cons1 .^ (1-sigma)) ./ (1 - sigma);
        ret1(cons1<0)=-Inf;
        ret1=reshape(ret1, [numel(z)*numel(a),1]);
        Q = makeQmatrix(pol_indx1, lzprob);
        for j=1:30
        vNew=ret1+beta*Q*V;
        v=vNew;
        j=j+1;
        end
        vfn=reshape(v,[numel(z) numel(a)]);
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:numel(z)
        E(:,:,i)=repmat(lzprob(i,:)*vfn, [numel(a) 1]);
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
   Mu = zeros(numel(z),numel(a));
  Mu(1, 4) = 1;
    %Mu=ones(nz, num_a)/(nz*num_a);
    % ITERATE OVER DISTRIBUTIONS
    mu_tol=1;
    while mu_tol >1e-7
          [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
          MuNew = zeros(size(Mu));
          for ii = 1:length(emp_ind)
              apr_ind = pol_indx(emp_ind(ii), a_ind(ii));
              MuNew(:, apr_ind) = MuNew(:, apr_ind) + (lzprob(emp_ind(ii),:)*Mu(emp_ind(ii),a_ind(ii)))';
          end
          mu_tol=max(max(abs(MuNew-Mu)));
          Mu=MuNew;
    end
    aggsav=sum(Mu*a');
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
% Euler equation error
c=bsxfun(@plus,bsxfun(@plus,-pol_fn,r*a),(z*w*l_bar)');
cf=c(:,pol_indx');
cf=reshape(cf,[numel(z) numel(a) numel(z)]);
for i=1:numel(z)
    c_prime(i,:)=lzprob(i,:)*cf(:,:,i);
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c_prime.^(-sigma)*r).*Mu));
