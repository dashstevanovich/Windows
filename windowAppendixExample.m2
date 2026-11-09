-- Verifying the example by computing the nonacyclic labels.
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

-- Verifying the example by first loading the nonacyclic labels:
-- Computing the nonacyclic labels is the most expensive part of this verification, so we saved
-- the images of the nonacyclic labels as strings in the folder exampleFiles so that one can load
-- these images. If one wishes to just check that these weights are too wide for L in the set Ls,
-- they may run the below code instead of the code above. We also saved the closed cells as
-- matrices consisting of their vertices, as this computation is also very expensive, however these
-- cells are only used to compute the nonacyclic images, and thus are not required in the below
-- computation.

load("windowAppendix.m2")
classMap = transpose matrix({{0,0,1},{0,0,1},{0,2,1},{0,2,1},{1,1,1},{2,1,1},{3,0,1}})
divisorVector = transpose matrix({{0,0,2,-5,7,5,-2}})
Iplus = {5,6}
Ls = {{},{0,1},{4},{2,3},{0,1,4},{2,3,4}}

weightsTooWideWrappedWithImages = L -> (
    imagesString = "exampleFiles/nonacyclicImages";
    scan(L, i -> imagesString = imagesString | toString(i));
    nonacyclicImages = value get imagesString;
    lambdaCone = getLambdaPlusMinusCone(classMap, Iplus, L);
    weightsTooWide(nonacyclicImages,lambdaCone,classMap, Iplus)
    );

apply(Ls, L -> weightsTooWideWrappedWithImages(L))

