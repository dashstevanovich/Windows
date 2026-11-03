load("windowAppendix.m2")
classMap = transpose matrix({{0,0,1},{0,0,1},{0,2,1},{0,2,1},{1,1,1},{2,1,1},{3,0,1}})
divisorVector = transpose matrix({{0,0,2,-5,7,5,-2}})
Iplus = {5,6}
Ls = {{},{0,1},{4},{2,3},{0,1,4},{2,3,4}}
divisorPolytope = getDivisorPolytope(classMap, divisorVector)
fList = apply(entries map(ZZ^(numcols classMap)), rowVec -> matrix({rowVec}))
allClosedCells = getClosedCells(divisorPolytope, fList)
weightsTooWideWrapped = L -> (
    nonacyclicLabels = getNonacyclicLabels(classMap, divisorPolytope, allClosedCells, L);
    nonacyclicImages = unique(apply(nonacyclicLabels, c -> classMap*c));
    lambdaCone = getLambdaPlusMinusCone(classMap, Iplus, L);
    weightsTooWide(nonacyclicImages,lambdaCone,classMap, Iplus)
    );

apply(L -> weightsTooWideWrapped(L))
