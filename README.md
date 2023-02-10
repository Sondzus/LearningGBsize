# LearningGBsize

This is code used for "Predicting the cardinality of a reduced Gröbner basis" [link coming as soon as preprint available] by Shahrzad Jamsidi, Eric Kang, and Sonja Petrović. 

The Jupyter notebooks [Vectors.ipynb](Vectors.ipynb) and [Features.ipynb](Features.ipynb) contain the code for learning and prediction described in the paper. 

The data files for 250,000 ideals are in the text files [zipped](4 text files data on 250000 binomial ideals.zip); the full data sets as described in the paper will be available shortly. 

Here is an examle on how to run the Macaulay2 code to generate samples. This is the same code we used in the paper.
Load the file [DataGeneratorScript.m2](DataGeneratorScript.m2) and then run the following commands in `M2`:

```
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

```
