# 1210

In this week , I implement a simulator for ACAM Tile , using LCA map to get the mismatch node , true leave are prediction of testing dataset , target leave are leaf node id of each path in a single tree.
using LCA_map[true_leave, target_leave] to index the mismatch feature. 


For the experiment in `code.r` , Try to using the two way , "naive" and "hamming distance" to choose the candidate path in a tile.
`naive` means sorting by path length
`hamming distance` means choosing the candidate path which have minimum hamming distance 


For the experiment in `code_hamming_revised.r` , I realize hamming distance is not a good upper bound estimation, so i try to using "union" to replace the "hamming distance"
`union` means try to find those path which have minimum '1' increasement in a bit vector , initial bit vector is using feature of the first path 

