-- This file contain starter code for taking a system of Polynomials in M2,
-- and auto generating the input files that multiregeneration.py can use
-- to solve it. 

--CURRENTLY SET UP FOR 6 BAR with variable pose/points

printingPrecision = 64

-- This is a helper function that adds commas in between
-- each pair of consecutive strings in a list. e.g. 
-- addCommas({"a","b","c"})
joinWithCommas = L -> (
    output = L;
    for i from 1 to length(L)-1 do 
        output = insert(2*i-1, ",", output);
    concatenate(output)
    )

-- STEP 1: Run whatever code you need to to get the system you want.
-- Here's an example system you might like to solve 
--R = CC[x,y,z]
--system = {x^2*y + y*z^2 + y*x + 2, y^2 + x^2 - z^2*y^2 + 1, x*y*z - 5}

--number of points
N = 2
--number of poses
M = 4 
--number of affine constraints we want to add
AFF_NUM = 0

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
    l12(d,G1,s,y1,t,z1)*l12(db,G1b,sb,y1b,tb,z1b))
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
angleconsts = flatten flatten for i from 1 to N list{
    {t_i*tb_i-1,s_i*sb_i-1}
}

angleconsts

--additional linear constraints
--takes in list of variables, outputs a linear constraint 
linGen = L ->(const = sum for i to #L-1 list(random(CC)*L#i);random(CC) + const)

--have to make angleVars a list
tList = flatten for i from 1 to N list(
    {t_i,tb_i,s_i,sb_i}
)

vargroup = join({G1,G1b,G2,G2b,G3,G3b,      --fixed points
    z1,z1b,z2,z2b,y1,y1b,y2,y2b},tList)  

--make affine constraint list
affineconsts = for i to AFF_NUM - 1 list(linGen(vargroup))

system = join(l1consts,l2consts,l3consts,angleconsts,affineconsts)



-- STEP 2: Choose variable groups.
-- Here are some example variable groups. Variable group 0 is {x,z} and variable group 1 is {y}
-- Recall that choosing different variable groups could speed up or slow down the computation.
-- What you want to look for when choosing variable groups is for variables in the same group (e.g. x and z)
-- not to be in the same terms. You can see that in the system above, x and z are rarely in the same term, so
-- the grouping below is a good guess for the best variable grouping.
--variableGroups = {{x,z}, {y}}


--CODE I ADDED
variableGroups = {vargroup}


-- If you don't want to worry about variable groups, then just make one
-- group contining all the variables by uncommenting the following line.

-- variableGroups = {{x,y,z}}

-- STEP 3: Write the file "bertiniInput_variables"

-- First we remove the file "bertiniInput_variables" if it exists already,
-- and start freash.
removeFile "bertiniInput_variables"

-- The following loop write each line of our new "bertiniInput_variables"
-- file one at a time
for group in variableGroups do (
    -- convert our list of variables into a list of strings
    variableStringList = for variable in group list toString(variable);
    
    -- Concatenate our list of strings together into one string,
    -- with commas in the middle.
    variableString = joinWithCommas(variableStringList);
    
    -- Add a line to the "bertiniInput_variables" file corresponding
    -- to the current variable group. See the M2 documentation for more info
    -- on how to use "<<"
    "bertiniInput_variables" << "variable_group " << variableString << ";"<<endl;
)
-- Now that we're done writing to "bertiniInput_variables", let's close it.
"bertiniInput_variables" << close

-- STEP 4: Write the file "bertiniInput_equations"

-- Similar to before, we remove the "bertiniInput_equations" file if 
-- it already exists.

removeFile "bertiniInput_equations"

-- Creat a list of function names f1,f2,..."
functionNames = for i from 1 to length(system) list "f"|toString(i);

-- Write the first line of the "bertiniInput_equations" file, which declares
-- the functions f1,f2,...
"bertiniInput_equations" << "function " << joinWithCommas(functionNames) << ";" <<endl

formattedSystem = for equation in system list replace("ii", "I", toString(equation))

-- Now we write a line defining each f0,f1,etc. to be the functions
-- in our system.
for i from 0 to length(system)-1 do
    "bertiniInput_equations" << functionNames_i|" = "|toString(formattedSystem_i)|";\n";

-- Close the file now that we are done
"bertiniInput_equations" << close


-- STEP 5: Write the file "inputFile.py"

-- Remove "inputFile.py" if it exists
removeFile "inputFile.py"


-- For each funtion, we want to list it's multidegree vector.
multidegreesOfEquationsInSystem = for function in system list 
(
    for group in variableGroups list
    (
       -- The degree of the equation "function" in the variable group
       -- "group" is given as the maximum of ...
       max(
           -- for each momomial in "function"...
       for monomial in terms(function) list 
       (
           -- the degree of the monomial in the variable group.
           -- The degree of a monomial in a group of variables is the sum
           -- of it's degree in each variable in the group.
           sum(group, variable -> degree(variable, monomial))
       )
       )
    )
)

-- The following code replaces "{" with "[" and "}" with "]". We need to do this
-- because M2 uses {} to define lists, whereas python uses [].
pythonFormattedMultidegrees = toString(multidegreesOfEquationsInSystem)
pythonFormattedMultidegrees = replace("\\{", "[",pythonFormattedMultidegrees)
pythonFormattedMultidegrees = replace("\\}", "]",pythonFormattedMultidegrees)

-- We now create the file "inputFile.py" as described in the tutorial.

-- This line adds the degree information
"inputFile.py" << "degrees = " << pythonFormattedMultidegrees << endl

"inputFile.py" << "verbose = 1" << endl

"inputFile.py" << "explorationOrder = \"depthFirst\"" << endl

"inputFile.py" << close
    
    
--Colin's fix!
removeFile "bertiniInput_trackingOptions"
"bertiniInput_trackingOptions" << "SECURITYLEVEL:1;" << endl << close
    
-- STEP 6: Running multiregeneration

-- Now that we've generated the input exit emacs or M2 and, run the command
-- "python2 ../multiregeneration.py"
-- from this directory.
 





