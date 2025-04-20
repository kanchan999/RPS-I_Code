% BOX37
% Compares SDMA, weighted sum, and NSGA-II for 2D polynomial MOOP.

function box37()
    rng(1); % For reproducibility

    %% (1) Solve Anchor Solutions with fmincon
    tic;
    anchorF1 = localMinimizeF(1); % Minimize f1
    x_star = anchorF1.x;
    f1_xstar = anchorF1.fvals(1);
    f2_xstar = anchorF1.fvals(2);

    anchorF2 = localMinimizeF(2); % Minimize f2
    y_star = anchorF2.x;
    f1_ystar = anchorF2.fvals(1);
    f2_ystar = anchorF2.fvals(2);
    anchorTime = toc;

    fprintf('Anchor Solutions:\n');
    fprintf('  x*: (%.4f, %.4f) => f1=%.4f, f2=%.4f\n', x_star(1), x_star(2), f1_xstar, f2_xstar);
    fprintf('  y*: (%.4f, %.4f) => f1=%.4f, f2=%.4f\n', y_star(1), y_star(2), f1_ystar, f2_ystar);

    % Nadir point for hypervolume
    nadir = [max(f1_xstar, f1_ystar), max(f2_xstar, f2_ystar)] * 1.1;

    %% (2) SDMA Implementation
    tic;
    [sdmaSolutions, sdmaFuncEvals] = run_sdma();
    sdmaTime = toc;
    sdmaObj = sdmaSolutions(:, 3:4);
    fprintf('SDMA completed in %.2f seconds\n', sdmaTime);

    %% (3) Weighted Sum Scalarization
    tic;
    numWeights = 100;
    weights = linspace(0, 1, numWeights);
    wsSolutions = [];
    wsFuncEvals = 0;
    for i = 1:numWeights
        lambda = weights(i);
        result = weighted_sum_scalarization(lambda);
        wsFuncEvals = wsFuncEvals + result.funcEvals;
        if ~isempty(result.x)
            fvals = evaluate_objectives_no_penalty(result.x);
            wsSolutions = [wsSolutions; [result.x, fvals]];
        end
    end
    if ~isempty(wsSolutions)
        [wsFronts, ~] = fast_non_dominated_sort(wsSolutions(:,3:4));
        wsSolutions = wsSolutions(wsFronts{1}, :);
        wsObj = wsSolutions(:, 3:4);
    else
        wsObj = [];
    end
    wsTime = toc;
    fprintf('Weighted Sum completed in %.2f seconds\n', wsTime);

    %% (4) Standalone NSGA-II
    tic;
    popSize = 100;
    maxGen = 100;
    initPop = rand(popSize, 2) * 5;
    nsga2Solutions = nsga2_standalone(initPop, maxGen);
    nsga2Obj = evaluate_objectives_no_penalty(nsga2Solutions);
    [nsga2Fronts, ~] = fast_non_dominated_sort(nsga2Obj);
    nsga2Solutions = nsga2Solutions(nsga2Fronts{1}, :);
    nsga2Obj = nsga2Obj(nsga2Fronts{1}, :);
    nsga2FuncEvals = popSize * maxGen * 3;
    nsga2Time = toc;
    fprintf('Standalone NSGA-II completed in %.2f seconds\n', nsga2Time);

    %% (5) Compute Metrics
    hv_sdma = compute_hypervolume(sdmaObj, nadir);
    hv_ws = compute_hypervolume(wsObj, nadir);
    hv_nsga2 = compute_hypervolume(nsga2Obj, nadir);

    spread_sdma = compute_spread(sdmaObj);
    spread_ws = compute_spread(wsObj);
    spread_nsga2 = compute_spread(nsga2Obj);

    [perc_sdma_dominated_by_ws, perc_ws_dominated_by_sdma] = compute_domination_percentage(sdmaObj, wsObj);
    [perc_sdma_dominated_by_nsga2, perc_nsga2_dominated_by_sdma] = compute_domination_percentage(sdmaObj, nsga2Obj);
    [perc_ws_dominated_by_nsga2, perc_nsga2_dominated_by_ws] = compute_domination_percentage(wsObj, nsga2Obj);

    %% (6) Create Comparison Tables
    comparisonTable = table(...
        {'SDMA'; 'Weighted Sum'; 'Standalone NSGA-II'}, ...
        [hv_sdma; hv_ws; hv_nsga2], ...
        [spread_sdma; spread_ws; spread_nsga2], ...
        [size(sdmaObj, 1); size(wsObj, 1); size(nsga2Obj, 1)], ...
        'VariableNames', {'Method', 'Hypervolume', 'Spread', 'NumSolutions'});
    disp('Comparison Table:');
    disp(comparisonTable);

    dominationTable = table(...
        {'SDMA vs WS'; 'SDMA vs NSGA-II'; 'WS vs NSGA-II'}, ...
        [perc_ws_dominated_by_sdma; perc_nsga2_dominated_by_sdma; perc_nsga2_dominated_by_ws], ...
        [perc_sdma_dominated_by_ws; perc_sdma_dominated_by_nsga2; perc_ws_dominated_by_nsga2], ...
        'VariableNames', {'Comparison', 'Perc_Dominated_by_First', 'Perc_Dominated_by_Second'});
    disp('Domination Percentage Table:');
    disp(dominationTable);

    performanceTable = table(...
        {'SDMA'; 'Weighted Sum'; 'Standalone NSGA-II'}, ...
        [sdmaTime; wsTime; nsga2Time], ...
        [sdmaFuncEvals; wsFuncEvals; nsga2FuncEvals], ...
        'VariableNames', {'Method', 'Time_s', 'FunctionEvaluations'});
    disp('Performance Table:');
    disp(performanceTable);

    %% (7) Visualize Pareto Fronts
    figure;
    hold on;
    grid on;
    if ~isempty(sdmaObj)
        scatter(sdmaObj(:, 1), sdmaObj(:, 2), 40, 'b', 'filled', 'DisplayName', 'SDMA');
    end
    if ~isempty(wsObj)
        scatter(wsObj(:, 1), wsObj(:, 2), 40, 'r', '^', 'DisplayName', 'Weighted Sum');
    end
    if ~isempty(nsga2Obj)
        scatter(nsga2Obj(:, 1), nsga2Obj(:, 2), 40, 'g', 's', 'DisplayName', 'Standalone NSGA-II');
    end
    xlabel('$$f_1(x) = \frac{(x_1 + x_2 - 7.5)^2}{4} + (x_2 - x_1 + 3)^2$$', ...
        'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$$f_2(x) = 0.4 \times (x_1 - 1)^2 + 0.4 \times (x_2 - 4)^2$$', ...
        'Interpreter', 'latex', 'FontSize', 16);
    title('Pareto Front Comparison', 'FontSize', 16);
    legend('Location', 'best', 'FontSize', 12);
    ax = gca;
    ax.LineWidth = 2.5;
    ax.FontSize = 16;
    hold off;
    saveas(gcf, 'pareto_front_comparison.png');
end

%% ================== SDMA Implementation ==================
function [globalNDSolutions, funcEvals] = run_sdma()
    rng(1);
    anchorF1 = localMinimizeF(1);
    anchorF2 = localMinimizeF(2);
    
    f1_lo = min(anchorF1.fvals(1), anchorF2.fvals(1));
    f1_hi = max(anchorF1.fvals(1), anchorF2.fvals(1));
    f2_lo = min(anchorF1.fvals(2), anchorF2.fvals(2));
    f2_hi = max(anchorF1.fvals(2), anchorF2.fvals(2));
    initialBox = [f1_lo, f2_lo, f1_hi, f2_hi];

    maxBoxIter = 6;
    refinedBoxes = refine_boxes_simple(initialBox, maxBoxIter);
    numBoxes = size(refinedBoxes, 1);
    funcEvals = 0;

    finalPopBoxes = cell(numBoxes, 1);
    for i = 1:numBoxes
        box = refinedBoxes(i, :);
        x0 = [];
        for tries = 1:5
            f1_rand = box(1) + rand*(box(3)-box(1));
            f2_rand = box(2) + rand*(box(4)-box(2));
            [x0_try, evals] = findPreimage_fmincon(f1_rand, f2_rand);
            funcEvals = funcEvals + evals;
            if ~isempty(x0_try)
                x0 = x0_try;
                break;
            end
        end

        if isempty(x0)
            initPop = rand(50, 2) * 5;
        else
            delta = 0.5;
            lb = max([0,0], x0 - delta);
            ub = min([5,5], x0 + delta);
            initPop = sampleNeighborhood(x0, lb, ub, 50);
        end

        [finalPop, evals] = nsga2_box_with_initpop(initPop, box, 50);
        funcEvals = funcEvals + evals;
        finalPopBoxes{i} = finalPop;
    end

    localSolutionsInBox = cell(numBoxes, 1);
    for i = 1:numBoxes
        box = refinedBoxes(i, :);
        pop = finalPopBoxes{i};
        if isempty(pop), continue; end
        fvals = evaluate_objectives_no_penalty(pop);
        inBoxMask = (fvals(:,1) >= box(1)) & (fvals(:,1) <= box(3)) & ...
                    (fvals(:,2) >= box(2)) & (fvals(:,2) <= box(4));
        subPop = pop(inBoxMask, :);
        subObj = fvals(inBoxMask, :);
        if ~isempty(subPop)
            [frontsLocal, ~] = fast_non_dominated_sort(subObj);
            localSolutionsInBox{i} = [subPop(frontsLocal{1}, :), subObj(frontsLocal{1}, :)];
        end
    end

    allSolutions = vertcat(localSolutionsInBox{:});
    if isempty(allSolutions), globalNDSolutions = []; return; end

    objAll = allSolutions(:, 3:4);
    [frontsGlobal, ~] = fast_non_dominated_sort(objAll);
    globalNDSolutions = allSolutions(frontsGlobal{1}, :);
end

%% ================== Helper Functions ==================
function [c, ceq] = nonlconG(z)
    x1 = z(1); x2 = z(2);
    g1 = -((x1 - 2)^3)/2 - x2 + 2.5;
    g2 = -x1 - x2 + 8*(x2 - x1 + 0.65)^2 + 3.85;
    c = [-g1; -g2];
    ceq = [];
end

function result = localMinimizeF(whichObj)
    bestCost = inf;
    bestX = [NaN, NaN];
    lb = [0, 0]; ub = [5, 5];
    for trial = 1:3
        x0 = lb + rand(1, 2) .* (ub - lb);
        opts = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 2000);
        switch whichObj
            case 1, fun = @objFunF1;
            case 2, fun = @objFunF2;
        end
        [xSol, costVal] = fmincon(fun, x0, [], [], [], [], lb, ub, @nonlconG, opts);
        if costVal < bestCost
            bestCost = costVal;
            bestX = xSol;
        end
    end
    result.x = bestX;
    result.fvals = [objFunF1(bestX), objFunF2(bestX)];
end

function val = objFunF1(z)
    x1 = z(1); x2 = z(2);
    val = ((x1 + x2 - 7.5)^2)/4 + (x2 - x1 + 3)^2;
end

function val = objFunF2(z)
    x1 = z(1); x2 = z(2);
    val = 0.4*(x1 - 1)^2 + 0.4*(x2 - 4)^2;
end

function fvals = evaluate_objectives_no_penalty(pop)
    N = size(pop, 1);
    fvals = zeros(N, 2);
    for i = 1:N
        x1 = pop(i,1); x2 = pop(i,2);
        fvals(i,1) = ((x1 + x2 - 7.5)^2)/4 + (x2 - x1 + 3)^2;
        fvals(i,2) = 0.4*(x1 - 1)^2 + 0.4*(x2 - 4)^2;
    end
end

function result = weighted_sum_scalarization(lambda)
    bestCost = inf;
    bestX = [NaN, NaN];
    lb = [0, 0]; ub = [5, 5];
    funcEvals = 0;
    for trial = 1:3
        x0 = lb + rand(1, 2) .* (ub - lb);
        opts = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 2000);
        fun = @(z) lambda*objFunF1(z) + (1-lambda)*objFunF2(z);
        [xSol, costVal, ~, output] = fmincon(fun, x0, [], [], [], [], lb, ub, @nonlconG, opts);
        funcEvals = funcEvals + output.funcCount;
        if costVal < bestCost
            bestCost = costVal;
            bestX = xSol;
        end
    end
    result.x = bestX;
    result.funcEvals = funcEvals;
end

function finalPop = nsga2_standalone(initPop, maxGen)
    pop = initPop;
    popSize = size(pop, 1);
    for gen = 1:maxGen
        obj = evaluate_population_no_box(pop);
        [fronts, rank] = fast_non_dominated_sort(obj);
        crowding = zeros(popSize, 1);
        for i = 1:length(fronts)
            crowding(fronts{i}) = calculate_crowding_distance(obj(fronts{i}, :));
        end
        matingPool = tournament_selection(pop, rank, crowding, popSize);
        offspring = [];
        lb = [0, 0]; ub = [5, 5];
        for i = 1:2:popSize
            p1 = matingPool(i, :);
            p2 = matingPool(mod(i,popSize)+1, :);
            [c1, c2] = sbx_crossover(p1, p2, lb, ub, 20);
            c1 = polynomial_mutation(c1, lb, ub, 20, 0.1);
            c2 = polynomial_mutation(c2, lb, ub, 20, 0.1);
            offspring = [offspring; c1; c2];
        end
        combinedPop = [pop; offspring];
        combinedObj = evaluate_population_no_box(combinedPop);
        [combinedFronts, ~] = fast_non_dominated_sort(combinedObj);
        newPop = [];
        curFront = 1;
        while size(newPop,1) + length(combinedFronts{curFront}) <= popSize
            newPop = [newPop; combinedPop(combinedFronts{curFront}, :)];
            curFront = curFront + 1;
        end
        if size(newPop,1) < popSize && curFront <= length(combinedFronts)
            remaining = popSize - size(newPop,1);
            lastFront = combinedFronts{curFront};
            distances = calculate_crowding_distance(combinedObj(lastFront, :));
            [~, sortIdx] = sort(distances, 'descend');
            selected = lastFront(sortIdx(1:remaining));
            newPop = [newPop; combinedPop(selected, :)];
        end
        pop = newPop;
    end
    finalPop = pop;
end

function obj = evaluate_population_no_box(pop)
    penaltyVal = 1e6;
    N = size(pop, 1);
    obj = zeros(N, 2);
    for i = 1:N
        x1 = pop(i,1); x2 = pop(i,2);
        f1_ = ((x1 + x2 - 7.5)^2)/4 + (x2 - x1 + 3)^2;
        f2_ = 0.4*(x1 - 1)^2 + 0.4*(x2 - 4)^2;
        g1 = -((x1 - 2)^3)/2 - x2 + 2.5;
        g2 = -x1 - x2 + 8*(x2 - x1 + 0.65)^2 + 3.85;
        if (x1<0 || x1>5 || x2<0 || x2>5 || g1<0 || g2<0)
            obj(i,:) = [penaltyVal, penaltyVal];
        else
            obj(i,:) = [f1_, f2_];
        end
    end
end

function refinedBoxes = refine_boxes_simple(initialBox, maxIter)
    rectangles_to_process = initialBox;
    retained_rectangles = [];
    for iteration = 1:maxIter
        new_rectangles_to_process = [];
        for r_idx = 1:size(rectangles_to_process,1)
            rect = rectangles_to_process(r_idx,:);
            f1_min = rect(1);  
            f2_min = rect(2);
            f1_max = rect(3);  
            f2_max = rect(4);
            f1_mid = (f1_min + f1_max)/2;
            f2_mid = (f2_min + f2_max)/2;
            subrects = [
                f1_min, f2_min, f1_mid, f2_mid;
                f1_min, f2_mid, f1_mid, f2_max;
                f1_mid, f2_min, f1_max, f2_mid;
                f1_mid, f2_mid, f1_max, f2_max
            ];
            for s_idx = 1:4
                sub = subrects(s_idx,:);
                corners = [
                    sub(1), sub(2);
                    sub(1), sub(4);
                    sub(3), sub(2);
                    sub(3), sub(4)
                ];
                cornerFeas = false(4,1);
                for c_idx = 1:4
                    cornerFeas(c_idx) = cornerFeasible_fmincon(corners(c_idx,1), corners(c_idx,2));
                end
                if ~( all(cornerFeas) || all(~cornerFeas) )
                    if iteration==maxIter
                        retained_rectangles = [retained_rectangles; sub];
                    else
                        new_rectangles_to_process = [new_rectangles_to_process; sub];
                    end
                end
            end
        end
        rectangles_to_process = new_rectangles_to_process;
        if isempty(rectangles_to_process), break; end
    end
    refinedBoxes = retained_rectangles;
end

function feasFlag = cornerFeasible_fmincon(f1_val, f2_val)
    nRestarts = 2;
    bestDist = inf;
    for attempt=1:nRestarts
        x0 = rand(1,2)*5;
        opts = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',500);
        [~, distVal] = fmincon(@(z)cornerDist(z,f1_val,f2_val), x0, [],[],[],[], [0,0],[5,5], @nonlconG, opts);
        if distVal<bestDist, bestDist=distVal; end
        if bestDist<1e-6, break; end
    end
    feasFlag=(bestDist<1e-6);
end

function val = cornerDist(z,f1t,f2t)
    x1=z(1); x2=z(2);
    f1_=((x1 + x2 -7.5)^2)/4 + (x2 - x1 +3)^2;
    f2_=0.4*(x1 -1)^2 +0.4*(x2 -4)^2;
    val=(f1_-f1t)^2 + (f2_-f2t)^2;
end

function [x0, funcEvals] = findPreimage_fmincon(f1_target, f2_target)
    bestDist = inf;
    bestSol = [];
    funcEvals = 0;
    for attempt = 1:2
        xInit = rand(1, 2) * 5;
        opts = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',500);
        [solX, distVal, ~, output] = fmincon(@(z)cornerDist(z,f1_target,f2_target), xInit, [],[],[],[], [0,0],[5,5], @nonlconG, opts);
        funcEvals = funcEvals + output.funcCount;
        if distVal < bestDist
            bestDist = distVal;
            bestSol = solX;
        end
    end
    if bestDist < 1e-6
        x0 = bestSol;
    else
        x0 = [];
    end
end

function [finalPop, funcEvals] = nsga2_box_with_initpop(initPop, box, maxGen)
    pop = initPop;
    popSize = size(pop,1);
    funcEvals = 0;
    for gen = 1:maxGen
        obj = evaluate_population_with_box(pop, box);
        [fronts, rank] = fast_non_dominated_sort(obj);
        crowding = zeros(popSize,1);
        for i = 1:length(fronts)
            crowding(fronts{i}) = calculate_crowding_distance(obj(fronts{i},:));
        end
        matingPool = tournament_selection(pop, rank, crowding, popSize);
        offspring = [];
        for i = 1:2:popSize
            p1 = matingPool(i,:);
            p2 = matingPool(mod(i,popSize)+1,:);
            [c1,c2] = sbx_crossover(p1,p2, [0,0], [5,5], 20);
            c1 = polynomial_mutation(c1, [0,0], [5,5], 20, 0.1);
            c2 = polynomial_mutation(c2, [0,0], [5,5], 20, 0.1);
            offspring = [offspring; c1; c2];
        end
        combinedPop = [pop; offspring];
        combinedObj = evaluate_population_with_box(combinedPop, box);
        funcEvals = funcEvals + size(combinedPop,1);
        [combinedFronts, ~] = fast_non_dominated_sort(combinedObj);
        newPop = [];
        curFront = 1;
        while size(newPop,1) + length(combinedFronts{curFront}) <= popSize
            newPop = [newPop; combinedPop(combinedFronts{curFront},:)];
            curFront = curFront + 1;
        end
        if size(newPop,1) < popSize
            remaining = popSize - size(newPop,1);
            lastFront = combinedFronts{curFront};
            distances = calculate_crowding_distance(combinedObj(lastFront,:));
            [~, sortIdx] = sort(distances, 'descend');
            newPop = [newPop; combinedPop(lastFront(sortIdx(1:remaining)),:)];
        end
        pop = newPop;
    end
    finalPop = pop;
end

function obj = evaluate_population_with_box(pop, box)
    penaltyVal = 1e6;
    N = size(pop,1);
    obj = zeros(N,2);
    for i = 1:N
        x1 = pop(i,1); x2 = pop(i,2);
        f1_ = ((x1 + x2 -7.5)^2)/4 + (x2 - x1 +3)^2;
        f2_ = 0.4*(x1 -1)^2 +0.4*(x2 -4)^2;
        g1 = -((x1 -2)^3)/2 - x2 +2.5;
        g2 = -x1 - x2 +8*(x2 - x1 +0.65)^2 +3.85;
        if (x1<0 || x1>5 || x2<0 || x2>5 || g1<0 || g2<0 || ...
            f1_<box(1) || f1_>box(3) || f2_<box(2) || f2_>box(4))
            obj(i,:) = [penaltyVal, penaltyVal];
        else
            obj(i,:) = [f1_, f2_];
        end
    end
end

function pop = sampleNeighborhood(x0, lb, ub, popSize)
    pop = lb + rand(popSize,2).*(ub - lb);
end

function hv = compute_hypervolume(obj, refPoint)
    if isempty(obj), hv = 0; return; end
    [~, idx] = sort(obj(:,1));
    obj_sorted = obj(idx,:);
    hv = 0;
    for i = 1:size(obj_sorted,1)-1
        hv = hv + (obj_sorted(i+1,1)-obj_sorted(i,1))*(refPoint(2)-obj_sorted(i,2));
    end
    hv = hv + (refPoint(1)-obj_sorted(end,1))*(refPoint(2)-obj_sorted(end,2));
end

function spread = compute_spread(obj)
    if isempty(obj) || size(obj,1) < 2, spread = NaN; return; end
    [~, idx] = sort(obj(:,1));
    obj_sorted = obj(idx,:);
    distances = sqrt(sum(diff(obj_sorted).^2,2));
    mean_dist = mean(distances);
    spread = sum(abs(distances - mean_dist)) / (size(obj_sorted,1)-1) / mean_dist;
end

function [perc_a, perc_b] = compute_domination_percentage(a, b)
    if isempty(a) || isempty(b), perc_a = NaN; perc_b = NaN; return; end
    count_a = 0;
    for i = 1:size(a,1)
        for j = 1:size(b,1)
            if dominates(b(j,:), a(i,:))
                count_a = count_a + 1;
                break;
            end
        end
    end
    count_b = 0;
    for i = 1:size(b,1)
        for j = 1:size(a,1)
            if dominates(a(j,:), b(i,:))
                count_b = count_b + 1;
                break;
            end
        end
    end
    perc_a = (count_a / size(a,1)) * 100;
    perc_b = (count_b / size(b,1)) * 100;
end

function [fronts, rank] = fast_non_dominated_sort(obj)
    N = size(obj,1);
    S = cell(N,1); n = zeros(N,1); rank = zeros(N,1); fronts = {[]};
    for p = 1:N
        S{p} = [];
        for q = 1:N
            if dominates(obj(p,:), obj(q,:))
                S{p} = [S{p}, q];
            elseif dominates(obj(q,:), obj(p,:))
                n(p) = n(p) + 1;
            end
        end
        if n(p) == 0
            rank(p) = 1;
            fronts{1} = [fronts{1}, p];
        end
    end
    i = 1;
    while ~isempty(fronts{i})
        Q = [];
        for p = fronts{i}
            for q = S{p}
                n(q) = n(q) - 1;
                if n(q) == 0
                    rank(q) = i + 1;
                    Q = [Q, q];
                end
            end
        end
        i = i + 1;
        fronts{i} = Q;
    end
    fronts = fronts(1:i-1);
end

function flag = dominates(a, b)
    flag = all(a <= b) && any(a < b);
end

function distance = calculate_crowding_distance(obj)
    [N, M] = size(obj);
    distance = zeros(N,1);
    for m = 1:M
        [vals, idx] = sort(obj(:,m));
        distance(idx(1)) = Inf;
        distance(idx(end)) = Inf;
        for i = 2:N-1
            distance(idx(i)) = distance(idx(i)) + (vals(i+1)-vals(i-1))/(max(vals)-min(vals));
        end
    end
end

function matingPool = tournament_selection(pop, rank, crowding, poolSize)
    N = size(pop,1);
    matingPool = zeros(poolSize, size(pop,2));
    for i = 1:poolSize
        a = randi(N); b = randi(N);
        if rank(a) < rank(b) || (rank(a) == rank(b) && crowding(a) > crowding(b))
            winner = a;
        else
            winner = b;
        end
        matingPool(i,:) = pop(winner,:);
    end
end

function [child1, child2] = sbx_crossover(p1, p2, lb, ub, eta)
    child1 = p1;
    child2 = p2;
    for i = 1:length(p1)
        if rand <= 0.5
            x1 = min(p1(i), p2(i));
            x2 = max(p1(i), p2(i));
            beta = 1 + 2*(x1 - lb(i))/(x2 - x1);
            alpha = 2 - beta^-(eta+1);
            if rand <= 1/alpha
                betaq = (rand*alpha)^(1/(eta+1));
            else
                betaq = (1/(2 - rand*alpha))^(1/(eta+1));
            end
            c1 = 0.5*((x1 + x2) - betaq*(x2 - x1));
            c2 = 0.5*((x1 + x2) + betaq*(x2 - x1));
            child1(i) = min(max(c1, lb(i)), ub(i));
            child2(i) = min(max(c2, lb(i)), ub(i));
        end
    end
end

function mutated = polynomial_mutation(child, lb, ub, eta, pm)
    mutated = child;
    for i = 1:length(child)
        if rand < pm
            delta1 = (mutated(i) - lb(i))/(ub(i) - lb(i));
            delta2 = (ub(i) - mutated(i))/(ub(i) - lb(i));
            r = rand;
            if r <= 0.5
                deltaq = (2*r + (1-2*r)*(1-delta1)^(eta+1))^(1/(eta+1)) - 1;
            else
                deltaq = 1 - (2*(1-r) + 2*(r-0.5)*(1-delta2)^(eta+1))^(1/(eta+1));
            end
            mutated(i) = mutated(i) + deltaq*(ub(i) - lb(i));
            mutated(i) = min(max(mutated(i), lb(i)), ub(i));
        end
    end
end