function out = amcmc(mcmc_def)
%AMCMC
%
% Adaptive MCMC
%
% Tiangang Cui, 17/Jan/2014

[def, out]      = process_input(mcmc_def);
v_curr          = def.init;
sigma           = def.sigma;

% initial covariance
switch def.cov_def
    case {'Fresh'}
        sample_cross= v_curr(:)*v_curr(:)';
        sample_sum  = v_curr;
        sample_num  = 1;
        M           = zeros(def.np,1);
        C           = eye(def.np);
        scale       = exp(sigma);
        L           = eye(def.np)*scale;
    case {'Restart'}
        M           = def.mean;
        C           = def.cov;
        sample_num  = def.sample_num;
        sample_sum  = M*sample_num;
        sample_cross= C*sample_num + sample_sum(:)*sample_sum(:)'/sample_num;
        scale       = exp(sigma)/sqrt(sum(C(:).^2));
        L           = chol(C+eye(def.np)*1E-3)'*scale;
end

switch def.proposal
    case {'MALA'}
        [mlpt_curr, mllkd_curr, mg_curr] = def.minus_log_post(v_curr);
    case {'RW'}
        [mlpt_curr, mllkd_curr] = def.minus_log_post(v_curr);
        mg_curr = 0;
end

% start MCMC
acc             = 0;
batch           = 0;

for i = 1:def.nstep
    switch def.proposal
        case {'MALA'}
            drift_curr  = -(C*mg_curr)*(scale^2/2);
            r_curr      = randn(def.np,1);
            v_next      = v_curr + drift_curr + L*r_curr;
            [mlpt_next, mllkd_next, mg_next] = def.minus_log_post(v_next);
            drift_next  = -(C*mg_next)*(scale^2/2); % drift_next is not re-used, because the adaptation
            r_next      = -( L\(drift_curr + drift_next) + r_curr );
            log_n2c     = - 0.5 * r_next(:)'*r_next(:);
            log_c2n     = - 0.5 * r_curr(:)'*r_curr(:);
            alpha       = (mlpt_curr - mlpt_next) + (log_n2c - log_c2n);
        case {'RW'}
            v_next      = v_curr + L*randn(def.np,1);
            [mlpt_next, mllkd_next] = def.minus_log_post(v_next);
            mg_next     = 0;
            alpha       = mlpt_curr - mlpt_next;
    end
    if log(rand) < alpha
        v_curr          = v_next;
        mg_curr         = mg_next;
        mlpt_curr       = mlpt_next;
        mllkd_curr      = mllkd_next;
        acc             = acc+1;
    end
    
    batch       = batch + 1;
    % update covariance
    sample_cross        = sample_cross + v_curr(:)*v_curr(:)';
    sample_sum          = sample_sum + v_curr;
    sample_num          = sample_num + 1;
    
    if  batch   == def.nbatch
        delta   = min(0.1,sqrt(def.nbatch/i));
        if (acc/def.nbatch) < def.rate
            sigma       = sigma - delta;
        else
            sigma       = sigma + delta;
        end
        batch   = 0;
        acc     = 0;
        % factorize new covariance
        M       = sample_sum/sample_num;
        C       = sample_cross/sample_num - M(:)*M(:)';
        scale   = exp(sigma)/sqrt(max(diag(C)));
        L       = chol(C)'*scale;
        
        out.k               = out.k+1;
        out.sigma(out.k,:)  = sigma;
    end
    
    % save
    out.grad2           = out.grad2 + mg_curr*mg_curr';
    out.v_samples(i,:)  = v_curr;
    out.lpt(i)          = mlpt_curr;
    out.llkd(i)         = mllkd_curr;
    out.mh(i)           = alpha;
    out.M               = M;
    out.C               = C;
end

out.grad2 = out.grad2/def.nstep;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [def, out] = process_input(def)


if ~isfield(def, 'proposal')
    def.proposal    = 'RW';
end
user.high_d         = false;
if strcmp(def.proposal, 'MALA')
    user.rate       = 0.58;
else
    user.rate       = 0.23;
end
[def, out]          = mcmc_input(def, user);

if ~isfield(def, 'cov_def')
    def.cov_def     = 'Fresh';
end

if strcmp(def.cov_def, 'Restart')
    if ~isfield(def, 'mean') || ~isfield(def, 'cov') || ~isfield(def, 'sample_num')
        disp('Error: no initial covariance defined');
    end
end

def.adapt_flag      = true;

out.grad2 = 0;

end
