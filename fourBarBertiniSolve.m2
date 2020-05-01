--This file solves a generic 4 bar linkage problem for 
--M poses and N positions. 
-- You can also scales the initial points, since it seems to make a difference numerically, but it should not conceptually.
-- You solve it through Bertini with your choice of precisions
--This files saves a short 1 line summey, the system that was solved and the solutions that were found in 3 different files 
--which can be found in teh same folder where it was run

--THis has been adjusted to homogenize the problem
restart;
printingPrecision =64
allowableThreads = 8;--allows it use multithreading, I thimlk we can increase it up to 8


needsPackage "Bertini";
-- I will homogenize with respect to G's Z's and Angles. 
--storeBM2Files =  temporaryFileName()
--mkdir storeBM2Files
fourBarBertiniSolve= (M,N,SCALE,METHOD) -> (
--generate an angle variable for each point
--(M,N,SCALE)=(1,4,1)
R = CC[G1,G1b,G2,G2b,hg]**CC[z1,z1b,z2,z2b,hz]**CC[Angle_1..Angle_N,Angleb_1..Angleb_N,rho]**CC[tau];
--generating random points
--function which gives pm 1
pm = i->(-1)^(random(2));


pointGen = i -> (SCALE/sqrt 2)*(pm()*random(RR)+ii*pm()*random(RR));


--generates our points
D = apply(1..N+M,pointGen);
DB = D/conjugate;


--define the angles
theta =toSequence apply(M,i->2*pi*random(RR));
--makes rotation for t in isotropic coordinates
toIsotropic = t-> cos(t)+ii*sin(t);


T = theta/toIsotropic ;--constants
TB = T/conjugate;


--Set up constraints
Z = (z1,z2);
ZB =(z1b,z2b);
G  =(G1,G2);
GB =(G1b,G2b);
A = Angle_1..Angle_N;
AB = Angleb_1..Angleb_N;
L = (d,theta,z,g) -> d+theta*z-g;
--print netList{Z,ZB,G,GB,A,AB};


--also a function of coordinates we want to use
createLij = (i,j)->(
     --this is a local variable so we use :=
     --appends variable angles to the fixed angles
     TA := T|A;
     TAB:= TB|AB;
     d:= D_i;
     db:=DB_i;
     z = Z_j;
     zb = ZB_j;
     gv = G_j;
     gvb = GB_j;
     theta = TA_i;
     thetab = TAB_i;
     complicatedL = L(d,theta,z,gv)*L(db,thetab,zb,gvb);
     if not member(theta,toList A) then (--print 1;
	 simpleL =rho*complicatedL) else(
--	 print 2;
	 coefAAB=diff(thetab,diff(theta,complicatedL));
--	 print coefAAB;
	 simpleL = complicatedL-(theta*thetab-1)*coefAAB;
  	 (affA,affGZ) = coefficients(simpleL,Variables=>A|AB);
	 --print (affA,affGZ);
	 simpleL = ((affA*(diagonalMatrix{1,1,rho})*affGZ))_(0,0)	;    	
	 );
     return simpleL
);
				
allL = apply(2,j->apply(M+N,i->createLij(i,j)));
--netList flatten allL

--arm j, point/pose i, list of all the length constraints
P = (i,j,ALL)->(
--because allL is list of two lists, each list for an arm
   ALL_j_0 - ALL_j_i
);
						
--constrains lnroms to be the same, arm j, point/pose i
createPij = apply((M+N-1),i->apply(2,j->P(i+1,j,allL)));
--netList createPij;

						
--Constraint that Angles should be conjugates
cutOff=0;
--angleConsts = flatten for i from 1 to N list if i>cutOff then (Angle_i -    rho) else Angleb_i-rho;

if METHOD===0 or METHOD ==="NAIVE" then 
angleConsts = flatten for i from 1 to N list  
(tau)*(random(CC)*Angle_i -  random CC*rho)*(random CC*Angleb_i-random CC*rho)+(1-tau)*(Angle_i*Angleb_i-1) else 
if METHOD ===1 or METHOD ==="MONIC" then
angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  random CC*rho)*(Angleb_i-random CC*rho)+(1-tau)*(Angle_i*Angleb_i-1) else 
if METHOD ===2 or METHOD ==="SYM" then(
    bfb:=for i from 1 to N list random CC;
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  bfb_(i-1)*rho)*(Angleb_i-bfb_(i-1)*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===3 or METHOD ==="SUPERSYM" then(
    b:= random CC;
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  b*rho)*(Angleb_i-b*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD === 4 or METHOD ==="CONJ" then(
    bfb=for i from 1 to N list random CC;
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  bfb_(i-1)*rho)*(Angleb_i-conjugate(bfb_(i-1))*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD === 5 or METHOD ==="SUPERCONJ" then(
    b= random CC;
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  b*rho)*(Angleb_i-conjugate(b)*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD === 6 or METHOD ==="MODULUS" then(
    bfalpha:=for i from 1  to N list random (RR);
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  exp(pi*ii*bfalpha_(i-1))*rho)*(Angleb_i- exp(-pi*ii*bfalpha_(i-1))*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===7 or METHOD ==="UNITPMONE" then(
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  1*rho)*(Angleb_i+1*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===8 or METHOD ==="UNITONE" then(
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  1*rho)*(Angleb_i-1*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===9 or METHOD ==="UNITii" then(
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  1.0*ii*rho)*(Angleb_i-1.0*ii*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===10 or METHOD ==="UNITPMii" then(
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  1.0*ii*rho)*(Angleb_i+1.0*ii*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    ) else 
if METHOD ===11 or METHOD ==="Origin" then(
    angleConsts = flatten for i from 1 to N list  (tau)*(Angle_i -  0*rho)*(Angleb_i+0*rho)+(1-tau)*(Angle_i*Angleb_i-1)
    );

					     
angleConsts = flatten angleConsts;		
--Linear constraints
affineConsts = flatten for i from 1 to (10-2*M-N) list{
    random(RR)*(G1+G1b)*hz+
    -random(RR)*ii*(G1-G1b)*hz+
    random(RR)*(G2+G2b)*hz+
    -random(RR)*ii*(G2-G2b)*hz+    
    random(RR)*(z1+z1b)*hg+
    -random(RR)*(z1-z1b)*hg+
    random(RR)*(z2+z2b)*hg+
    -random(RR)*(z2-z2b)*hg+
    hz*hg*random RR    
};	
		
pList = flatten createPij;
--length should be (M+N-1)*2
#pList == (M+N-1)*2;

homPoly = (f,H)->(
    D := degree f;
    newPoly := 0;
    sum apply(terms f,m->m*product apply(#H,i->(H_i)^(D_i-(degree m)_i)))
    );
ourSystem = join(affineConsts,pList,angleConsts);
ourSystem=ourSystem/(f->homPoly(f,{hg,hz,rho}));

makeB'InputFile(storeBM2Files,
    HomVariableGroup=>{{G1,G1b,G2,G2b,hg},{z1,z1b,z2,z2b,hz},    
    	toList(Angle_1..Angle_N)|toList(Angleb_1..Angleb_N)|{rho}},
    BertiniInputConfiguration=>{UseRegeneration=>1,MPType=>2,SecurityLevel=>1},
    B'Constants=>{tau=>1},
    B'Polynomials=>ourSystem
    );
runBertini(storeBM2Files);
unfilteredPreCircle = importSolutionsFile(storeBM2Files,NameSolutionsFile=>"nonsingular_solutions");
print("Number unfilteredPreCircle sols",#unfilteredPreCircle);
--We say a solution is degenerate if G1=G2, G1b=G2b,z1=z2,z1=z1b;
--We filter out the degenerate solutions by a tolerance
solsPreCircle ={};
tolerance=1e-8;
scan(unfilteredPreCircle,p->if {p_0-p_2,p_1-p_3,p_4-p_6,p_5-p_7}/abs//min>1e-10 then solsPreCircle=solsPreCircle|{p} );
print("Number unfilteredPreCircle sols",#unfilteredPreCircle);
print("Number solsPreCircle sols",#solsPreCircle);

moveB'File(storeBM2Files,"input","input1");
moveB'File(storeBM2Files,"bertini_session.log","session1");
moveB'File(storeBM2Files,"nonsingular_solutions","ns_PreCircle_Unfiltered_"|M|"_"|N|METHOD);
writeStartFile(storeBM2Files,solsPreCircle);
writeParameterFile(storeBM2Files,{1},NameParameterFile=>"start_parameters");
writeParameterFile(storeBM2Files,{0},NameParameterFile=>"final_parameters");
makeB'InputFile(storeBM2Files,
    HomVariableGroup=>{{G1,G1b,G2,G2b,hg},{z1,z1b,z2,z2b,hz},    
    	toList(Angle_1..Angle_N)|toList(Angleb_1..Angleb_N)|{rho}},
    ParameterGroup=>{tau},
    BertiniInputConfiguration=>{MPType=>2,SecurityLevel=>1,ParameterHomotopy=>2},
    B'Polynomials=>ourSystem
    );
runBertini(storeBM2Files);

unfilteredPostCircle = importSolutionsFile(storeBM2Files,NameSolutionsFile=>"nonsingular_solutions");
--We filter out the degenerate solutions by a tolerance as before. 
solsPostCircle ={};
tolerance=1e-8;
scan(unfilteredPostCircle,p->if {p_0-p_2,p_1-p_3,p_4-p_6,p_5-p_7}/abs//min>1e-10 then solsPostCircle=solsPostCircle|{p}) ;
print("Number unfilteredPostCircle sols",#unfilteredPostCircle);
print("Number solsPostCircle sols",#solsPostCircle);
moveB'File(storeBM2Files,"nonsingular_solutions","ns_PostCircle_Unfiltered_"|M|"_"|N|METHOD);


toString (M,N,SCALE,"Bertini System") << ourSystem <<close;
--toString (M,N,SCALE,"Bertini Solutions") << solsPostCircle <<close;
toString(M,N)<<toString ("Case: ", M, N, SCALE, METHOD,
     "Number of unfilteredPreCircle: ", #unfilteredPreCircle,
     "Number of solsPreCircle: ", #solsPreCircle,
     "Number of unfilteredPostCircle: ", #unfilteredPostCircle,
     "Number of solsPostCircle: ", #solsPostCircle)<<close;

("Case: ", M, N, SCALE, "Number of solutions using Bertini: ", #solsPostCircle) 
);




print(fourBarBertiniSolve(4,1,1,"MODULUS"));
print("a")


for m in (2,1) do (
    for n from 0 to 10-2*m do (
    print (fourBarBertiniSolve(m,n,1,"NAIVE"));
    ));













print "Gold"
print storeBM2Files
storeBM2Files = temporaryFileName()
mkdir storeBM2Files
print(fourBarBertiniSolve(4,1,1,"NAIVE"));--

print(fourBarBertiniSolve(4,1,1,"NAIVE"));--
"""
(Number unfilteredPreCircle sols, 60)
(Number PreCircle sols, 50)
(Number solsPostCircle sols, 48)
(Case: , 4, 1, 1, Number of solutions using Bertini: , 48)
"""

print(fourBarBertiniSolve(3,1,1,"NAIVE"));--
print(fourBarBertiniSolve(2,1,1,"NAIVE"));--
print(fourBarBertiniSolve(2,3,1,"NAIVE"));--
print(fourBarBertiniSolve(1,6,1,"NAIVE"));

print(fourBarBertiniSolve(1,3,1,"NAIVE"));-- 234 and two? degenerate from 254 precircle




