# 260520

In this week , I added four embedded dataset from UCI repo. 
`Energy` : https://archive.ics.uci.edu/dataset/242/energy+efficiency
`breast` : https://archive.ics.uci.edu/dataset/17/breast%2Bcancer%2Bwisconsin%2Bdiagnostic
`heart`  : https://archive.ics.uci.edu/dataset/45/heart%2Bdisease
`parkinsons` : https://archive.ics.uci.edu/dataset/189/parkinsons%2Btelemonitoring

and `code.r` evaluate the total performance of each datasets, result show in `experiment_result` , this file evaluate the read energy of a testing data and total write energy of each method on each dataset.

by using `analyze.r` summarize the `experiment_result` and calulate "total read energy"( read_energy_pj * testing_dataset_size ) and "total write energy" ( write_energy_pj ) of each method on dataset.

fake_total_energy_reduction_ratio = 0.1 means fake match saving 10% of total energy compare to rewrite.

`dataset_setup.r` preprocess all of the embedded datasets locate in ($data_dir)

 
