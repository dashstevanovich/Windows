loadPackage("SimplicialComplexes",Reload => true)

getDivisorPolytope = method()
getDivisorPolytope(Matrix, Matrix) := (classMap, divisorVector) -> (
    zeroVector := matrix(mutableMatrix(ZZ,numrows classMap,1));
    imageOfM := affinePreimage(classMap,convexHull(zeroVector));
    positiveOctant := polyhedronFromHData(-1*map(ZZ^(numcols classMap)), divisorVector);
    return intersection(imageOfM,positiveOctant)
    );

slicePolytopeByHyperplanes = method()
slicePolytopeByHyperplanes(Polyhedron,Matrix) := (P,f) -> (
    fImage := affineImage(f, P);
    intervalBounds := (entries vertices fImage)_0;
    minValue := min(apply(intervalBounds, vertex -> floor(vertex)));
    maxValue := max(apply(intervalBounds, vertex -> ceiling(vertex)));
    slicesInHyperplane := for i from minValue to maxValue list (
	hyperplane := polyhedronFromHData(f || -f, matrix({{i},{-i}}));
	sliceInHyperplane := intersection(P,hyperplane);
	if dim sliceInHyperplane == -1 then continue else sliceInHyperplane
	);
    slicesInBetween := for i from minValue to maxValue-1 list (
	region := polyhedronFromHData(f || -f, matrix({{i+1},{-i}}));
	intersection(P,region)
	);
    return join(slicesInHyperplane, slicesInBetween)
    );

collectUnique = method()
collectUnique(List) := (inList) -> (
    outList := {};
    scan(inList, C -> if not any(outList, Q -> Q == C) then outList = append(outList, C));
    return outList
    );

getClosedCells = method()
getClosedCells(Polyhedron,List) := (P, fList) -> (
    allCells = {P};
    for f in fList do (
	cellsByf := for Q in allCells list (
	    slicePolytopeByHyperplanes(Q,f)
	    );
	allCells = collectUnique(flatten(cellsByf));
	);
    return allCells
    );

getCellCeiling = method()
getCellCeiling(Polyhedron) := (P) -> (
    x := interiorPoint(P);
    xEntries := (entries transpose x)_0;
    xEntriesCeiling := apply(xEntries, xi -> ceiling(xi));
    return transpose matrix({xEntriesCeiling})
    );

getNerve = method()
getNerve(Matrix, List, Polyhedron) := (cellCeiling, L, divisorPolytope) -> (
    n := ambDim divisorPolytope;
    Lcomp := select(toList(0..n-1), i -> not isMember(i,L));
    identityMatrix := map(ZZ^n);
    cellCeilingMinusOne := cellCeiling - matrix(table(n,1,(i,j) -> 1));
    rectangularRegionInequalities := submatrix(identityMatrix,Lcomp,) || -1*submatrix(identityMatrix,Lcomp,);
    rectangularRegionValues :=  submatrix(cellCeiling,Lcomp,) || -1*submatrix(cellCeilingMinusOne,Lcomp,);
    SdeltaLClosed := intersection(polyhedronFromHData(rectangularRegionInequalities,rectangularRegionValues),divisorPolytope);
    subsetsOfLcomp := subsets(Lcomp);
    lowerFacetRegions := new HashTable from apply(subsetsOfLcomp,
	J -> J => polyhedronFromHData(submatrix(identityMatrix,J,) || -1*submatrix(identityMatrix,J,),
	    submatrix(cellCeilingMinusOne,J,) || -1*submatrix(cellCeilingMinusOne,J,)));
    return select(subsetsOfLcomp, J -> dim(intersection(SdeltaLClosed, lowerFacetRegions#J)) >= 0)
    );

getNonacyclicLabels = method()
getNonacyclicLabels(Matrix,Polyhedron,List,List) := (classMap,divisorPolytope,allClosedCells,L) -> (
    n := numcols classMap;
    S := ZZ[x_0..x_(n-1)];
    SVars := vars S;
    allCellCeilings := unique(apply(allClosedCells, Q -> getCellCeiling(Q)));
    allNerves := apply(allCellCeilings, cellCeiling -> getNerve(cellCeiling,L,divisorPolytope));
    allNerveMonomials := apply(allNerves,
	nerve ->  apply(nerve,
	    J -> promote(product(flatten entries(submatrix(SVars,J))),S)));
    allNerveComplexes := apply(allNerveMonomials, nerveMonomials -> simplicialComplex(nerveMonomials));
    nonacyclicLabelIndices := positions(allNerveComplexes, Delta -> HH Delta != 0);
    return apply(nonacyclicLabelIndices, i -> allCellCeilings_i)
    );

getLambdaPlusMinusCone = method()
getLambdaPlusMinusCone(Matrix, List, List) := (classMap, Iplus, Izero) -> (
    n := numcols classMap;
    Iminus := select(toList(0..n-1), i -> not isMember(i, Iplus) and not isMember(i, Izero));
    greaterThanIneqs := transpose (submatrix(classMap, Iplus) | submatrix(classMap, Izero));
    lessThanIneqs := transpose (submatrix(classMap, Iminus) | submatrix(classMap, Izero));
    return coneFromHData(greaterThanIneqs || -1*lessThanIneqs)
    );

weightsTooWide = method()
weightsTooWide(List,Cone,Matrix,List) := (nonacyclicLabels,lambdaCone,classMap,Iplus) -> (
    CMaxCones := apply(nonacyclicLabels, label1 -> intersection(lambdaCone,
	    coneFromHData(matrix(apply(nonacyclicLabels,
		    label2 -> flatten entries transpose (label1-label2))))));
    CMinCones := apply(nonacyclicLabels, label1 -> intersection(lambdaCone,
	   coneFromHData(matrix(apply(nonacyclicLabels,
		    label2 -> flatten entries transpose (label2-label1))))));
    betaIplus := submatrix(classMap,Iplus)*matrix(table(length(Iplus),1,(i,j) -> 1));
    labelIterator := 0..(length(nonacyclicLabels)-1);
    scan(labelIterator, i -> scan(labelIterator, j ->
	    if not contains(coneFromHData(
		    transpose(nonacyclicLabels_i - nonacyclicLabels_j-betaIplus)),
		intersection(CMaxCones_i,CMinCones_j)) then return false));
    return true
    );

end
restart
