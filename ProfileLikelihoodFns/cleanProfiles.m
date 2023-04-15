function out = cleanProfiles(out,threshold)

% this function will take in a collection of profiles (out) along with a
% threshold and "clean" them. Clean means make sure the value does not
% cross the threshold more than once in either direction AND trim the ends
% so that we only have the 95% CI and not however far over that we stepped.

dbstop if naninf % help with cleaning profiles manually

sz = size(out);
npars = size(out,1);
out = reshape(out,npars,[]);
n_abm_vecs = size(out,2);

for i = 1:n_abm_vecs
    for j = 1:npars
        try % if this returns without error, then just trim the profile
            out{j,i} = trimProfile(out{j,i},threshold);
        catch % but if there's an error, then there's multiple crosses, rerun the optimization to smooth it out
            [I1,I2,~] = findMaxCrosses(out{j,i}(end,:),threshold);
            [i1,~] = find(I1'-I2==1); % I1 contains indices where the value began to drop; if these come right after an increase over the threshold, then these are likely outliers we want to ignore
            if numel(i1)==1 % then there was one time a jump occurred, ignore this
                out{j,i}(:,I1(i1)) = [];
                out{j,i} = trimProfile(out{j,i},threshold);
            else
                fprintf("Unknown error encountered in trimming profile %d at parameter %d. \n   You can fix manually now and then continue.\n",i,j)
                1/0; % cause execution to pause here to inspect issue
            end
        end
    end
end

out = reshape(out,sz);

dbclear if naninf
