-- ******************************************************************
-- Data generator script that is copied all to one file for the purpose of sharing the code on GitHub.  
-- ******************************************************************

-* 
 Random binomial data is written to SEVERAL FILES. 
 The filename for this data set is of- the format "RandomBinomialDataSet.3vars.deg4.sampleSize10". 

 The format of the file is: 
 1st line contiains parameters for the data generated,
 2nd line is a list of lists of binomial exponents, printed as string from array,

 Other output files that go with each sample: 
 filename.features.txt (first line =  header line) 
 filename.gbSizes.txt
 filename.gbMaxDeg.txt 

 TimeOut option: if the gb computation times out for a particular ideal, then the corresponding line in the file gbsizes will be "0" and gbMaxDeg will be "-1". 

*- 
-- ******************************************************************
-- PREREQUISITES 
-- ******************************************************************
--loadPackage "RunExternalM2";  -- for timing out the gb runs 



-- ******************************************************************
--      randomBinomials
-- ****************************************************************** -- load"randomBinomialIdeals.m2"
-- make a list "mons" of all mons of degree d in ring R; 
-- then select two random indices, indPlus and indMinus,
-- and make a binomial out of mons_indPlus - mons_indMinus,
-- and add k such binomials to a list. 
-- 
-- ExcludeZeroBinomials=>true, controls whether zero binomials are allowed
-- if Homogeneous=>true then: 
-- OUTPUT: list of k homogeneous binomials (of format M1-M2) of degree d in ring R.
-- if Homogeneous=>false then: 
-- OUTPUT: list of k (non-homogeneous) binomials (of format M1-M2) of degree at most d in ring R.
---------------------------------------------------------------------------------------------
randomBinomials = method(TypicalValue => List, Options=>{Homogeneous=>false, ExcludeZeroBinomials=>true})
randomBinomials(PolynomialRing,ZZ,ZZ) := List => o -> (R,maxDegree,k) -> ( 
    if o.Homogeneous then mons = flatten entries basis(maxDegree,R) else (
	mons = {};
	scan(0..maxDegree, d-> mons = append(mons, flatten entries basis(d,R)));
	mons = flatten mons
	);
    mons = drop(mons,1); -- no 0.  TO DO: enable this later; it may produce monomials not true binomials.
    binomials = {};
    scan(k, i-> (
	    indPlus = random(0,#mons-1);
	    indMinus = random(0,#mons-1);
	    if (o.ExcludeZeroBinomials==true) then(
	    	while (indPlus==indMinus) do(
	    	    indPlus = random(0,#mons-1);
	      	    indMinus = random(0,#mons-1);
	      	    );
		);
	    --
	    binomials = append(binomials, mons_indPlus-mons_indMinus);
	    )
    	);
    flatten binomials
    )
-- but the list may include zeros. 
-- ******************************************************************


-- ******************************************************************
-- GENERATE GB training DATA 
-- run this section of code once to get BOTH
-- the random binomial ideals saved as a text file,
-- and the size of each minimal GB saved in a different text file. 
-- ******************************************************************
generateGBdata = method(TypicalValue => List, Options=>{TimeLimit=>null,Homogeneous=>false,InfoLine=>true,InfoLineCompact=>false,ExcludeZeroBinomials=>true,MonOrder=>null,SaveFeatures=>null})
generateGBdata(ZZ,ZZ,ZZ,ZZ) := List => o -> (numVars,maxDegree,binomialsInEachSample,sampleSize) -> ( 
    if o.MonOrder===null then  S = ZZ/32003[x_0..x_(numVars-1)] else  S=ZZ/32003[x_0..x_(numVars-1),MonomialOrder=>o.MonOrder];
--    filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|".sampleSize"|toString sampleSize|"."|toString binomialsInEachSample|"binomialsEach"|".txt";
--    filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|"."|toString binomialsInEachSample|"binomialsEach.txt";
    -*
always same filename! no longer adding timestamp!! 
    *-
    -- generate a random binomial set, write it to file f, 
    -- then compute gb write size gb and max deg of gb, 
    -- save any custom features as well  as the input  data hardcoded features, 
    -- and  then continue to next random set.
	scan(sampleSize,i-> (
		-- generate a new random binomial set: 
		bins = randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>o.Homogeneous,ExcludeZeroBinomials=>o.ExcludeZeroBinomials);
		assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
		-- save exponents to the open file f: 
		expos =  toString apply(sort bins,b->apply(exponents b,monexpo->monexpo));	 -- added `sort` to ensure unique rep of the data set
		-- but first get rid ouf outer {} bc they serve no purpose:
	    	--f = openOut concatenate(filename|".txt");
		f = openOutAppend concatenate(filename|".txt");
    		--f = openOutAppend filename;
		f << replace(" ","",substring(1,#expos-2,expos))<<endl;
       		close f;
		
		-- if there's a time limit we want to impose on each GB computation: 
    	    	if o.TimeLimit===null then (
		    I := ideal toList bins;
		    myGB := gens gb I;
		    (gbSize,gbMaxDeg) := (# flatten entries myGB, (max (degrees  myGB)_1)_0);
		    gbRemove (I); -- not to keep the huge GBs in memory - we no longer need them! 
		    I = symbol I;
		    myGB = symbol myGB;
		    ) else ( 
		    	childRun = runExternalM2("Gstvari.m2","Gstvari",(toString bins,S),PreRunScript=>"ulimit -t 2");
			-- childRun#value -- this is the output of the function Gstvari; childRun is a hashtable w/ bunch of statistics from the child M2 process..
			if childRun#value=!=null then (gbSize,gbMaxDeg) = childRun#value else  (gbSize,gbMaxDeg) = (0,-1); 
		    );

		gbf := openOutAppend concatenate(filename,".gbSizes.txt");
		gbf <<  gbSize << endl;
		close gbf;
		--ADDING FILE FOR SAVING MAX TOTAL DEGREE OF A GB ELEMENT:  
		gbdegf := openOutAppend concatenate(filename,".gbMaxDeg.txt");
		-- the  next  line computes the maximum of the list  of degrees  of the elements of the GB:  (_1 to handle whatever else M2 outputs) 
    	    	gbdegf << gbMaxDeg << endl;
		close gbdegf;		
		--end of adding. 
		
	    -- FEATURES:
	    -- The ideal way (pun not intended) is if all these features are in the same file and separated by commas.
	    -- To clarify: each sample is a separate line, and the features are separated by commas. 
		-- generator degree statistics (min, max, mean, variance): 
	    	degreesOfGenerators := bins/degree; --shorthand for "apply(bins,degree)" 
	    	mean := (sum degreesOfGenerators)/#degreesOfGenerators; 
		variance :=  sum apply(degreesOfGenerators, x-> (mean_0-x_0)^2)/#degreesOfGenerators;
		-- 1) built-in hard-coded features: 
		featurefile = openOutAppend concatenate(filename,".features.txt");
		featurefile <<  (min degreesOfGenerators)_0 <<","<< (max degreesOfGenerators)_0 <<","<< toString mean_0 <<"," << toString variance; 
    	    	-- 2) any custom (on-demand) features from optional input: 
		if o.SaveFeatures =!= null then (
		    -- assuming SaveFeatures = {function1,function2}  and each can be applied to List:
		    features  = apply(o.SaveFeatures,fn-> fn  bins); 
		    apply(length o.SaveFeatures,k-> (
			    featurefile << ","<< features_k; 
			    )
			);
		    );
    	    	featurefile << endl; -- done with this sample; 
		close featurefile; 
		
		bins = symbol bins;
		expos = symbol expos;
	    )
	);
	-- Eventually one may wish totime out the operations and if gb doesn't complete we can save a 0 or a -1 in the gb sizes file: 
	--  http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.14/share/doc/Macaulay2/RunExternalM2/html/_run__External__M2.html
	--  http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.14/share/doc/Macaulay2/RunExternalM2/html/_suggestions_spfor_spusing_sp__Run__External__M2.html
    print("Saved data to files starting with "|filename);
)


end  -- stop reading this file on Load. The code below has been copied to README on https://github.com/Sondzus/LearningGBsize 

-- ******************************************************************
-- This is the data generation code to run after the methods above have been loaded
-- ******************************************************************

-- preset parameters:
numVars = 5
maxDegree = 15
binomialsInEachSample = 5 

filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|"."|toString binomialsInEachSample|"binomialsEach.txt";
    
-- store the preset parameters for this data set on the first line of the data file: 
parameters = "numVars = "|toString numVars|", maxDegree = "|toString maxDegree|", binomialsInEachSample = "|toString binomialsInEachSample|", MonomialOrder = default";
	    f = openOut filename;
	    f << parameters << endl;
	    close f;
	    featurefile = openOut concatenate(filename,".features.txt");
	    featurefile << "min_deg_gen, max_deg_gen, mean_deg_gen, var_deg_gen, numgens, dim, deg, reg" << endl; --- Edit this code if you'd like to work with other features and write down what they are.
	    close featurefile;


generateGBdata(numVars,maxDegree,binomialsInEachSample,sampleSize=500, Homogeneous=>false,  
    SaveFeatures => {numcols@@mingens@@ideal,dim@@ideal,degree@@ideal,regularity@@ideal}  
--   ,  InfoLine=>false, InfoLineCompact=>true
    )
    -- ,    TimeLimit => 1) ---optional input for limiting the time we allow a GB computation to run for each sample point.

-- ******************************************************************
