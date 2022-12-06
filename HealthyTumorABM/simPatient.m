function M = simPatient(M)

M = finishParameterSetup_Patient(M);

%% initialize inputs
M.t = 0;
M.i = 1; % step index
M = initializeEnvironment(M);
M = initializeTracked(M);

%% initialize figure
if M.plot_pars.plot_fig
    M = initializeFigure(M);
end

%%
if M.save_pars.dt < Inf
    M = finishParameterSetup_Saves(M);
    M = saveInitialModelData(M);
else
    M.next_save_time = Inf;
end

if M.plot_pars.make_movie
    M = initializeMovie(M);
end

for ei = 1:M.events.n % event index
    M = finishParameterSetup_Event(M,M.events.times(ei)-M.t);
    M = simForward(M);

    extendTimeSeriesPlots(M,ei)

    switch M.events.event_index(ei)

        case 1 % start tumor
            M = initializeTumor(M);

        case Inf % end of simulation
            break;
    end
end

if M.save_pars.dt < Inf
    M = saveFinalModelData(M);
end

if M.plot_pars.make_movie
    close(M.vid)
end
