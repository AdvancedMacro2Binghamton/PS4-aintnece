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
        pol_indx_1=pol_indx-1;
        pol_indx_1(pol_indx_1==0)=1;
        pol_indx_2=pol_indx+1;
        pol_indx_2(pol_indx_2==501)=500;
        a_new=[a(pol_indx_1),a(pol_indx),a(pol_indx_2)];
        for i=1:numel(z)
            for j=1:numel(a)
              pol_indx_1(j,:,i)=pol_indx_1(j,:,i)+numel(a)*(j-1)+numel(a)^2*(i-1);
              pol_indx_2(j,:,i)=pol_indx_2(j,:,i)+numel(a)*(j-1)+numel(a)^2*(i-1);
              pol_indx_m(j,:,i)=pol_indx(j,:,i)+numel(a)*(j-1)+numel(a)^2*(i-1);
            end
        ret_t(:,:,i)=ret(:,:,i)';
        E_t(:,:,i)=E(:,:,i)';
        end
        ret_1=ret_t(pol_indx_1);
        E_1=E_t(pol_indx_1);
        ret_2=ret_t(pol_indx_2);
        E_2=E_t(pol_indx_2);
        ret_m=ret_t(pol_indx_m);
        E_m=E_t(pol_indx_m);
        ret_new=[ret_1,ret_m,ret_2];
        E_new=[E_1,E_m,E_2];
        b=[1,2,3];
        b1=interpn(b,5);
        i=1;
        for i=1:numel(z)
           ret_new1(:,:,i)=(interp1(b,ret_new(:,:,i)',b1))';
           E_new1(:,:,i)=(interp1(b,E_new(:,:,i)',b1))';
           a_new1(:,:,i)=(interp1(b,a_new(:,:,i)',b1))';
        end
        v_new=ret_new1+beta*E_new1;
        [vmax_new, pol_indx_new]=max(v_new, [], 2);
        vmax_new=permute(vmax_new, [3 1 2]);
        
        
        v_tol = max(max(abs(vmax_new-v_guess)));
        v_guess=vmax_new;
    end
     % KEEP DECSISION RULE
    pol_indx_new=permute(pol_indx_new, [3 1 2]);
        for i=1:numel(z)
            for j=1:numel(a)
              pol_indx_new1(i,j)=pol_indx_new(i,j)+65*(j-1)+500*65*(i-1);
            end
        a_new_t(:,:,i)=a_new1(:,:,i)';
        end
    
     pol_fn=a_new_t(pol_indx_new1);
     for i=1:numel(z)
     d(:,:,i) = abs(bsxfun(@minus, pol_fn(i,:),a'));
     [a_min, a_indx]=min(d(:,:,i) ,[],1);
     p(i,:)=a_indx;
     di(i,:)=pol_fn(i,:)-a(a_indx);
     end
   P= di./(a(2)-a(1))+p;
   
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
