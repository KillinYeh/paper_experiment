# Experiment
In this week , we want to find the upper bound of power saving , when we consider the similarity of paths
we can close the maximum bit line to saving power , and this perspective is not considerd in the relative paper.

So we made a experiment to compare naive vs proposed (considerd similarity ).



## Step
In the `code.r` I try to 
1. using those features which used by top 128 frequently feature to decide the first path in the tile.
2. find the paths which are head 127 shortest hamming distance within top 128 feature.
Repeat step 1 2 until all of the paths in the forest are included.


In the `code_hamming.r` , I try to 
1. Using all of the feature and find the path which have the smallest path length to be first path, ==instead of top 128 frequently used feature==
2. find the paths which are head 127 shortest hamming distance within ==all of the feature==.
3. calculate the number of four situations ( WL/BL , open/closed).  
