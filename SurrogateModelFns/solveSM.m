function sm_data = solveSM(sm,p,tt,C,D,condition_on_previous,resample_t)

persistent resample;
persistent sm_pre_processor;
persistent sm_solver;
persistent sm_post_processor;

if isempty(resample)
    resample = ~isempty(resample_t);
end

if isempty(sm_pre_processor)
    if isfield(sm,"pre_processor")
        sm_pre_processor = sm.pre_processor;
    else
        sm_pre_processor = @identityFunction;
    end
end

if isempty(sm_solver)
    switch sm.type
        case "ode"
            if condition_on_previous
                sm_solver = @solveODEWithExtend;
            else
                sm_solver = @solveODE;
            end
            if ~isfield(sm,"processor")
                if resample
                    sm_post_processor = @(raw_sm_output) deval(resample_t,raw_sm_output)';
                else
                    sm_post_processor = @(raw_sm_output) deval(tt,raw_sm_output)';
                end
            else
                sm_post_processor = sm.processor;
            end
    end

end

[sm,p,tt,C,~] = sm_pre_processor(sm,p,tt,C,D);
raw_sm_output = sm_solver(sm,p,tt,C);
sm_data = sm_post_processor(raw_sm_output);

end

function raw_sm_output = solveODE(sm,p,tt,C)
raw_sm_output = sm.solver(@(t,x) sm.fn(t,x,p,C),[sm.t0,tt(end)],sm.y0);
end

function raw_sm_output = solveODEWithExtend(sm,p,tt,C,D)
raw_sm_output = sm.solver(@(t,x) sm.fn(t,x,p,C),[sm.t0,tt(1)],sm.y0);
for i = 2:length(tt)
    raw_sm_output = odextend(raw_sm_output, [], tt(i), D.A(i-1,:));
end
end

function varargout = identityFunction(varargin)
    varargout = varargin;
end