#read the .irem files
wgas2 <- read.table("wgas2.irem", header=FALSE, stringsAsFactors=FALSE)
wgas3 <- read.table("wgas3_step2.irem", header=FALSE, stringsAsFactors=FALSE)

#get the IDs (in this case, first and second column are the same) 
wgas2_ids <- wgas2$V1
wgas3_ids <- wgas3$V1

#find IDs in wgas2 but NOT in wgas3 (excluded in approach 1, included in approach 2)
in_wgas2_not_wgas3 <- setdiff(wgas2_ids, wgas3_ids)

# print results
print(in_wgas2_not_wgas3)

