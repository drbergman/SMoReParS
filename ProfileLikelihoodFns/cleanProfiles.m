function profiles = cleanProfiles(profiles,threshold)

% this function will take in a collection of profiles (out) along with a
% threshold and "clean" them. Clean means make sure the value does not
% cross the threshold more than once in either direction AND trim the ends
% so that we only have the 95% CI and not however far over that we stepped.

sz = size(profiles);
npars = size(profiles,1);
profiles = reshape(profiles,npars,[]);
n_abm_vecs = size(profiles,2);

for i = 1:n_abm_vecs
    for j = 1:npars
        if i==178 && j == 4
            disp('')
        end
        profiles{j,i} = cleanThisProfile(profiles{j,i},i,j,threshold);
    end
end

profiles = reshape(profiles,sz);

end

function profile = cleanThisProfile(profile,i,j,threshold)

try % if this returns without error, then just trim the profile
    profile = trimProfile(profile,threshold);
catch % but if there's an error, then there's multiple crosses, remove extraneous times it rose above the threshold
    [inds_before_drop,inds_before_rise,max_val] = findMaxCrosses(profile(end,:),threshold);
    [i1,~] = find(inds_before_drop'-inds_before_rise==1); % I1 contains indices where the value began to drop; if these come right after an increase over the threshold, then these are likely outliers we want to ignore
    if ~isempty(i1) % then there was at least one time a jump occurred, ignore these
        profile(:,inds_before_drop(i1)) = [];
        profile = cleanThisProfile(profile,i,j,threshold);
    else
        % plot unclean profile
        f=figureOnRight;
        plot(profile(j,:),profile(end,:),"Marker",".")
        yline(max_val)

        [ST,current_ws_ind] = dbstack();
        stopText = strsplit(sprintf('in %s at %d', ST(current_ws_ind).name, ST(current_ws_ind).line + 3));
        dbstop(stopText{:});
        fprintf("Unknown error encountered in trimming profile %d at parameter %d. \n   You can fix manually now and then continue.\n",i,j)

        if false %#ok<*UNRCH> % suppress warning that the following lines cannot be reached, they are here to make them easily copy-pasted into Command Window
            % if just multiple one-step peaks over max_val (this can no longer be the fix since I have changed the check above from numel(i1)==1 to ~isempty(i1), meaning I will remove all one-step peaks)
            profile(:,inds_before_drop(i1)) = [];

            % if there's a single, multi-step peak over max_val
            profile(:,(inds_before_rise(i2_ind)+1):inds_before_drop(i1_ind)) = [];
            
        end
        hold on;
        plot(profile(j,:),profile(end,:))
        profile = cleanThisProfile(profile,i,j,threshold);
        close(f)
    end
end

end
