function MDProfile = performMultiDimProfile(files,objfn_constants,par_vals,force_serial)


sz = reshape(cellfun(@numel,par_vals),1,[]);
N = prod(sz);
npars = length(sz);
if ~isfield("objfn_constants","p_setup_fn")
    objfn_constants.p_setup_fn = @(p) p;
end

load(files.data_file,"t","D","C","cohort_size");
m = size(D,1);
D = reshape(D,m,[]); % string out all the cohorts along the 2nd dim
n_abm_vecs = size(D,2); % number of ABM parameter vectors used

if isfield(files,"previous_profile_file")
    load(files.previous_profile_file,"MDProfile")
    MDProfile = reshape(MDProfile,n_abm_vecs,N);
else
    MDProfile = zeros(n_abm_vecs,N);
end

field_sz = size(D(1).A);
n = prod(field_sz);
A = arrayify(D,"A",true);
A = reshape(A,[field_sz,m,n_abm_vecs]);
S = arrayify(D,"S",true);
S = reshape(S,[field_sz,m,n_abm_vecs]);

fn = objfn_constants.fn;
fn_opts = objfn_constants.fn_opts;
weights = reshape(objfn_constants.weights,1,[]);
weights = weights/sum(weights);
use_weights = numel(unique(weights))>1;

P = allCombos(par_vals{:},'matlab');

if force_serial
    for i = 1:N
        p = P(i,:)';
        sim_data = zeros([field_sz,m]);
        for j = 1:m
            sim_data(:,:,j) = fn(p,t,C{j},fn_opts);
        end
        if use_weights
            MDProfile(:,i) = -0.5*weights*(n*log(2*pi()) + reshape(sum(S.^2,1:2) + sum(((sim_data - A)./S).^2,1:2,"omitnan"),m,n_abm_vecs));
        else
            MDProfile(:,i) = -0.5*(n*log(2*pi()) + sum(S.^2,1:3) + sum(((sim_data - A)./S).^2,1:3,"omitnan"));
            if any(MDProfile(:,i)>0)
                error("positive ll")
            end
        end
    end
else
    if isempty(gcp('nocreate'))
        ppool = parpool("Threads"); %#ok<NASGU>
    else
        ppool = gcp; %#ok<NASGU>
    end
    C = parallel.pool.Constant(C);
    parfor i = 1:N
        p = P(i,:)';
        sim_data = zeros([field_sz,m]);
        for j = 1:m
            sim_data(:,:,j) = fn(p,t,C.Value{j},fn_opts); %#ok<PFBNS> (suppresses warning about fn bearing an overhead cost
        end
        if use_weights
            MDProfile(:,i) = -0.5*weights*(n*log(2*pi()) + reshape(sum(S.^2,1:2) + sum(((sim_data - A)./S).^2,1:2,"omitnan"),m,n_abm_vecs));
        else
            MDProfile(:,i) = -(0.5/m)*(n*log(2*pi()) + sum(S.^2,1:3) + sum(((sim_data - A)./S).^2,1:3,"omitnan")); % divide by m to have weights sum to 1 as I do when weights are not identical
        end
    end
end

MDProfile = reshape(MDProfile,[cohort_size,sz]);