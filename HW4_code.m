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
a_max=7;
num_a=500;

d=1;
while d>1e-5
% Discretization of assets

a=linspace(a_min,a_max,num_a);

% VFI
k_min=a_bar;
k_max=a_max;

dis=1;

while abs(dis)>=0.01
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

    % SET UP INITITAL DISTRIBUTION

    Mu=ones(size(pol_fn))/numel(pol_fn);
    % ITERATE OVER DISTRIBUTIONS
    
    MuNew = zeros(size(Mu));
    mu_tol=1;
    while mu_tol>1e-06
        [emp_ind, a_ind, mass] = find(Mu);
        MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (lzprob(emp_ind(ii), :) * mass(ii))';
    end
    mu_tol=max(max(abs(MuNew-Mu)));
    Mu=MuNew;
    end
    dis=k_guess-sum(sum(Mu.*pol_fn));
    if dis>0
        k_max=k_guess;
    else
        k_min=k_guess;
    end
end
d=Mu(end,end);
a_max=a_max+1;
end
% Policy functions graph

figure;
plot(a,pol_fn);
title('Policy Function with different Z');
legendcell=cellstr(num2str(z'));
legend(legendcell,'Location','northwest');

% Wealth Distribution

wd=bsxfun(@times,Mu,a);
figure;
plot(a,wd);
title('Wealth Distribution with Different Z');
legend(legendcell,'Location','northwest');

% Lorenz Curve & Gini Coefficient
w=repmat(reshape(a,[num_a 1]),numel(z),1);
mu=reshape(Mu,[numel(pol_fn) 1]);
c=cumsum(sortrows([mu mu.*w w],3));
L=bsxfun(@rdivide,c,c(end,:))*100;
wi=repmat(w,1,numel(pol_fn));
wj=repmat(w',numel(pol_fn),1);
g=sum(sum(abs(wi-wj)))/sum(w)/numel(pol_fn)/2;
figure;
area(L(:,1),L(:,2));
axis([0 100 0 100]);
axis('auto y')
hold on;
fplot(@(x) x,':');
title(['Gini Coefficient=',num2str(g)]);
