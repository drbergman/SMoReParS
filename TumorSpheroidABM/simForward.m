function M = simForward(M)

if M.Nsteps == 0 || isempty(M.tumor) % nothing to simulate here
    return;
end

%% prepare new tracked values
names = fieldnames(M.tracked);
for i = 1:length(names)
    M.tracked.(names{i}) = cat(1,M.tracked.(names{i}),...
        zeros([M.Nsteps,size(M.tracked.(names{i}),2:ndims(M.tracked.(names{i})))]));
end

%% iterations
for i = 1:M.Nsteps

    if M.NT > M.pars.max_tumor_size
        names = fieldnames(M.tracked);
        for ni = 1:length(names)
            colons = repmat({':'},1,ndims(M.tracked.(names{ni}))-1);
            M.tracked.(names{ni})(M.i+1:end,colons{:}) = NaN;
        end
        break
    end

    M.t = M.t + M.dt;
    M.i = M.i+1;

    % set these counts here to be updated as cells proliferate or are found
    % to be in interphase
    M.tracked.phase_cell_hours(M.i,:) = [0,M.NT * M.dt,0]; % [interphase,ready to proliferate,did proliferate]

    M = updateTumor(M);

    %% clean up tumor stuff
    M = removeApoptotic(M);

    %% advance tumor proliferation timers
    M.tumor(:,M.I.proliferation_timer) = max(0,M.tumor(:,M.I.proliferation_timer)-M.dt); % update time until tumor cells can proliferate for all that did not proliferate (1) or were just created (0) because those timers were updated in updateTumor

    %% update these vals
    M.NT = size(M.tumor,1);

    %% update tracked values
    M = updateTracked(M);

    %% plot
    if M.plot_pars.plot_fig && mod(M.i,M.plot_pars.plot_every)==M.plot_pars.plot_offset
        plotFunction_EndStep(M)
    end

    %% movie
    if M.plot_pars.plot_fig && M.plot_pars.make_movie
        writeVideo(M.vid,print(M.fig.handle,'-r100','-RGBImage'))
    end

    %% save any "big" data
    if M.t >= M.next_save_time - 0.5 * M.dt % then it appears that this is the closest time to the desired save time
        M = saveModelData(M);
    end

end %%end of for