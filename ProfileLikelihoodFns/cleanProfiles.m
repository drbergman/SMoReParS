function out = cleanProfiles(out,threshold)

% this function will take in a collection of profiles (out) along with a
% threshold and "clean" them. Clean means make sure the value does not
% cross the threshold more than once in either direction AND trim the ends
% so that we only have the 95% CI and not however far over that we stepped.

sz = size(out);
npars = size(out,1);
out = reshape(out,npars,[]);
n_abm_vecs = size(out,2);

% [ST,current_ws_ind] = dbstack();

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
                % plot unclean profile
                max_val = min(out{j,i}(end,:)) + threshold;
                f=figureOnRight;
                plot(out{j,i}(j,:),out{j,i}(end,:))
                yline(max_val)

                [ST,current_ws_ind] = dbstack();
                stopText = strsplit(sprintf('in %s at %d', ST(current_ws_ind).name, ST(current_ws_ind).line + 3));
                dbstop(stopText{:});
                fprintf("Unknown error encountered in trimming profile %d at parameter %d. \n   You can fix manually now and then continue.\n",i,j)

                if false %#ok<*UNRCH> % suppress warning that the following lines cannot be reached, they are here to make them easily copy-pasted into Command Window
                    % if just multiple one-step peaks over max_val
                    out{j,i}(:,I1(i1)) = [];

                    % if there's a single, multi-step peak over max_val
                    out{j,i}(:,(I2(i2_ind)+1):I1(i1_ind)) = [];

                    % check on cleaned profile
                    out{j,i} = trimProfile(out{j,i},threshold);
                    hold on;
                    plot(out{j,i}(j,:),out{j,i}(end,:))
                end
                
                close(f)
            end
        end
    end
end

out = reshape(out,sz);
