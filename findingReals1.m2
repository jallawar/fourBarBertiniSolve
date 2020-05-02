--this file is trying to find reals for 4 bar linkage with some specified points
restart

--needsPackage "PHCpack"
needsPackage "Bertini"

--functions which takes number of points and outputs some information
partitionSolve = N->(

R = CC[G1x,G1y,Px_0..Px_(N-1), Py_0..Py_(N-1)]; 

--boolean for which package to use
BERT = true;
--boolean for if you want to keep running even if errors are thrown
TRAINING_WHEELS = true;


--generating x,y points
dx = toSequence apply(N, i->(-1)^(random(2))*random(RR)*1); 
dy = toSequence apply(N, i->(-1)^(random(2))*random(RR)*1);


--this constrains that x_i is the same distance from the fitted point and G1
topCircleLens = apply(N, i->(
	(Px_i - dx_i)^2+(Py_i - dy_i)^2 
));

--constrains lengths to be the same
topCircleConsts = apply(N-1,i->(
    topCircleLens#(N-1)-topCircleLens#i	
));

bottomCircleLens = apply(N, i->(    
	(Px_i - G1x)^2+(Py_i - G1y)^2 
));

--constrains lengths to be the same
bottomCircleConsts = apply(N-1,i->(
    bottomCircleLens#(N-1)-bottomCircleLens#i	
));


--for affine constraints
--vargrps = join({{G1x},{G1y}},apply(N,i->{Px_i}),apply(N,i->{Py_i}))
vargrps = join({{G1x,G1y}},apply(N,i->{Px_i,Py_i}));
--generates random real number
ranReal = i->(-1)^(random(2))*random(RR);


--generates K affine constraints
--affineGen = K-> apply(K,i->ranReal()+sum apply(varList, j->ranReal()*j))

--this sums variables by their little group and gives them the same coefficients
affineGen = (K,varGroups)-> apply(K,i->ranReal() + sum apply(varGroups,j->ranReal()*sum j));


firstSystem = join(topCircleConsts,bottomCircleConsts,affineGen(4,vargrps));


-*
this might be the worst thing ever to write in a program, but for some reason
PHCpack shits itself if you run it once, but if you keep running it it gets it 
eventually. don't ask me why  
*-

if TRAINING_WHEELS then(
    hasntRun = true;
    while hasntRun do(
        try(if BERT then sols = bertiniZeroDimSolve firstSystem else 
	    sols = solveSystem firstSystem) then(hasntRun = false) else(
	    hasntRun = true)
	)
) else sols = solveSystem firstSystem;



sols = for i in sols list(coordinates i);
-- reals check

tol = -10; --larger negative number means stricter tolerance
realCheck = i-> (abs imaginaryPart i)<10^(tol);


realSols = for i in sols list (
    reals = for j in i list(
        if realCheck(j) then true else break;
    );
    if #reals < #i then continue else i        
);


--code for other half 

-*
So we set up a similar thing as the first half, but now with a 
both circle constraints and a dot product constraint between z1 and z2
*-

--sets number of affine constraints for second part
if N < 4 then AFF_NUM = 2 else if N == 4 then AFF_NUM = 1 else AFF_NUM = 0;


--function which gets real solutions to the other half..hopefully!
--function of index for solution in realSols
otherHalf = k ->( 
    
    sol = realSols#k;

    R2 = CC[G2x,G2y,Qx_0..Qx_(N-1),Qy_0..Qy_(N-1)];
   
    --this constrains that x_i is the same distance from the fitted point and G1
    topCircleLens = apply(N, i->(
	(Qx_i - dx_i)^2+(Qy_i - dy_i)^2 
    ));

    --constrains lengths to be the same
    topCircleConsts = apply(N-1,i->(
        topCircleLens#(N-1)-topCircleLens#i	
    ));

    bottomCircleLens = apply(N, i->(    
	(Qx_i - G2x)^2+(Qy_i - G2y)^2 
    ));

    --constrains lengths to be the same
    bottomCircleConsts = apply(N-1,i->(
        bottomCircleLens#(N-1)-bottomCircleLens#i	
    ));

    --dot product constraint between z1 and z2
    z1x = apply(N, i->sol_(2+i) - dx_i);
    z1y = apply(N, i->sol_(2+N+i) - dy_i);
    
    dotprods = apply(N, i->z1x_i*(Qx_i - dx_i)+z1y_i*(Qy_i-dy_i));
    
    --fixes all dot products to be same
    dotProdConsts = apply(N-1,i-> dotprods#(N-1) -dotprods#i);
    
    --for affine constraints
    --varList2 = join({G2x,G2y},apply(N,i->Qx_i),apply(N,i->Qy_i));
    varList2 = join({{G2x,G2y}},apply(N,i->{Qx_i,Qy_i}));
    --generates K affine constraints
    --affineGen = K-> apply(K,i->(-1)^(random(2))*random(RR)*50 + sum apply(varList2, j->random(RR)*(-1)^(random(2))*j));
    
    secondSystem = join(topCircleConsts,bottomCircleConsts,dotProdConsts,affineGen(AFF_NUM,varList2));
    
    if TRAINING_WHEELS then(
        hasntRun = true;
        while hasntRun do(
            try(if BERT then sols2 = bertiniZeroDimSolve secondSystem else
	        sols2 = solveSystem secondSystem) then( hasntRun = false) else(
	        hasntRun = true)    
        )
    ) else sols2 = solveSystem secondSystem;
    
    sols2 = for i in sols2 list(coordinates i);
    
    realSols2 = for i in sols2 list (
     	reals = for j in i list(
            if realCheck(j) then true else break;
        );
        if #reals < #i then continue else i        
    );

    if #realSols2 == 0 then print("no real solutions found to other half") else print("gotcha!");
    return realSols2
);


secondHalfReals = for i to #realSols - 1 list(
    otherHalf(i)
); 

--# of solutions for each real first half
--solsPerSol = for i in secondHalfReals list(#i);

 --sum solsPerSol
 
 list(realSols,secondHalfReals)

)


--stop skipping!



start3 = cpuTime()
elapsedTime partitionSolve(3)
end3=cpuTime()
print end3 - start3

start4 = cpuTime()
elapsedTime partitionSolve(4)
end4=cpuTime()
end4-start4



--rounds = for i to 9 list(partitionSolve(3))





