# 0301

## "MCP"
Adding the new way to find out the "similar" path , we treat path as bit vectors, and we try to find the most common prefix by radix sort (O(N*F)).
## radix + union
MCP only deal with prefix , its will let the postfix inpredictable and cause bit line can't not be closed .
So to the postfix , I decided to choose the candidate by union in small window size (setup 512 here).

  
