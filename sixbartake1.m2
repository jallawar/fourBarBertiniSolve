--
--This file will try to solve a Stephenson 3 six bar linkage which has 3 fixed points G and two triangles
--we're constraining every angle at every point
--

restart
needsPackage "PHCpack"


--number of points
N = 0
--number of poses
M = 5 
--number of affine constraints we want to add
AFF_NUM = 2

R' = CC[G1,G1b,G2,G2b,G3,G3b,      --fixed points
    z1,z1b,z2,z2b,y1,y1b,y2,y2b,   --triangle arms
    t_1..t_N,s_1..s_N, --angle variables for points   
    tb_1..tb_N,sb_1..sb_N]

--try poses that angles for both 
--generates a random complex point 
pointGen = i-> 100*(random(RR)*(-1)^(random(2))+random(RR)*(-1)^(random(2))*ii)
--takes an angle, outputs cos(t)+isin(t)
angleGen = t->cos(t)+ii*sin(t)

--generates our points
D = toSequence apply(M+N,pointGen)
DB = D/conjugate

--generates our angles for S, T
theta = toSequence apply(M,i->random(RR)*2*pi) 
T = theta/angleGen
TB = T/conjugate
--appends variables
T= T|(t_1..t_N)
TB = TB|(tb_1..tb_N)

theta2 = toSequence apply(M,i->random(RR)*2*pi)
S = theta2/angleGen
SB = S/conjugate
S = S|(s_1..s_N)
SB = SB|(sb_1..sb_N)

--functions which makes l constraints
l12 = (d,g,s,y,t,z)->d-g+s*y+t*z -- works for l1 and l2 just different y
l3 = (d,g,t,z)->d-g+t*z

--ith point, L1
createL1 = i->(
        d:= D_i;
	db:= DB_i;
	t := T_i;
	tb := TB_i;
	s := S_i;
	sb := SB_i;
	l12(d,G1,s,y1,t,z1)*l12(db,G1b,sb,y1b,tb,z1b)
)
--ith point, L2
createL2 = i->(
        d:= D_i;
	db:= DB_i;
	t := T_i;
	tb := TB_i;
	s := S_i;
	sb := SB_i;
	l12(d,G2,s,y2,t,z1)*l12(db,G2b,sb,y2b,tb,z1b)
)
--ith point, L1
createL3 = i->(
        d:= D_i;
	db:= DB_i;
	t := T_i;
	tb := TB_i;
	l3(d,G3,t,z2)*l3(db,G3b,tb,z2b)
)

--lists of lengths 
l1s = apply(M+N,createL1)
l2s = apply(M+N,createL2)
l3s = apply(M+N,createL3)
--function to make sure all lengths are the same
lengthConstrainer = (i,LIST)-> LIST_0 - LIST_i

--makes length constraints
l1consts = apply(M+N-1,i->lengthConstrainer(i+1,l1s))
l2consts = apply(M+N-1,i->lengthConstrainer(i+1,l2s))
l3consts = apply(M+N-1,i->lengthConstrainer(i+1,l3s))

--constraints for angle variables
angleconsts = flatten for i from 1 to N list{
    {t_i*tb_i-1,s_i*sb_i-1}
}

--additional linear constraints
--takes in list of variables, outputs a linear constraint 
linGen = L ->(
    const = sum for i to #L-1 list(random(CC)*L#i);
    random(CC) + const
)

vargroup = {G1,G1b,G2,G2b,G3,G3b,      --fixed points
    z1,z1b,z2,z2b,y1,y1b,y2,y2b}   --triangle arms

--make affine constraint list
affineconsts = for i to AFF_NUM - 1 list(
    linGen(vargroup)    
)

ourSystem = join(l1consts,l2consts,l3consts,angleconsts,affineconsts)

sols = solveSystem(ourSystem)




