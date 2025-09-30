# 251001
Replace average "row" utilization by calculate average "tile" utilization
tile utilization means consider fixed row and columns in one tile , here set rows = columns = 256.


## Issue
If columns << tile_size , its will waste so many space ( if dataset have 8 feature totally and tile_size = 256 , in this case , waste column = 256 - 8)
so if allow , change column_size to 2^(floor(lg(Nfeat))) , try to minimize the waste space.
fully code see in code2.r .
