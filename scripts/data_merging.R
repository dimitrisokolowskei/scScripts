# Brief Description ----
# Here I import and create an object Seurat for each individual sample. Finally, I merged each one into the same Seurat object and metadata.
# Importing Libraries ----
library(SeuratObject)
library(Seurat)

# Data importing
# 2047 samples ----
F47G0 <- Read10X(data.dir = "/work/dimitri/single/2047_matrix/filtered_2047_gex0") # Filtered 2047 Gex0
F47G0_obj <- CreateSeuratObject(counts = F47G0, min.features = 100)

F47G7 <- Read10X(data.dir = "/work/dimitri/single/2047_matrix/filtered_2047_gex7")
F47G7_obj <- CreateSeuratObject(counts = F47G7, min.features = 100) 

F47G21 <- Read10X(data.dir = "/work/dimitri/single/2047_matrix/filtered_2047_gex21")
F47G21_obj <- CreateSeuratObject(counts = F47G21, min.features = 100) 

F47G28 <- Read10X(data.dir = "/work/dimitri/single/2047_matrix/filtered_2047_gex28")
F47G28_obj <- CreateSeuratObject(counts = F47G28, min.features = 100)  

# 2049 samples ----
F49G0 <- Read10X(data.dir = "/work/dimitri/single/2049_matrix/filtered_2049_gex0") 
F49G0_obj <- CreateSeuratObject(counts = F49G0, min.features = 100)

F49G7 <- Read10X(data.dir = "/work/dimitri/single/2049_matrix/filtered_2049_gex7") 
F49G7_obj <- CreateSeuratObject(counts = F49G7, min.features = 100)

F49G21 <- Read10X(data.dir = "/work/dimitri/single/2049_matrix/filtered_2049_gex21") 
F49G21_obj <- CreateSeuratObject(counts = F49G21, min.features = 100)

F49G28 <- Read10X(data.dir = "/work/dimitri/single/2049_matrix/filtered_2049_gex28") 
F49G28_obj <- CreateSeuratObject(counts = F49G28, min.features = 100)

# 2051 samples ----
F51G0 <- Read10X(data.dir = "/work/dimitri/single/2051_matrix/filtered_2051_gex0") 
F51G0_obj <- CreateSeuratObject(counts = F51G0, min.features = 100)

F51G7 <- Read10X(data.dir = "/work/dimitri/single/2051_matrix/filtered_2051_gex7") 
F51G7_obj <- CreateSeuratObject(counts = F51G7, min.features = 100)

F51G21 <- Read10X(data.dir = "/work/dimitri/single/2051_matrix/filtered_2051_gex21") 
F51G21_obj <- CreateSeuratObject(counts = F51G21, min.features = 100)

F51G28 <- Read10X(data.dir = "/work/dimitri/single/2051_matrix/filtered_2051_gex28") 
F51G28_obj <- CreateSeuratObject(counts = F51G28, min.features = 100)

# 2052 samples ----
F52G0 <- Read10X(data.dir = "/work/dimitri/single/2052_matrix/filtered_2052_gex0") 
F52G0_obj <- CreateSeuratObject(counts = F52G0, min.features = 100)

F52G7 <- Read10X(data.dir = "/work/dimitri/single/2052_matrix/filtered_2052_gex7") 
F52G7_obj <- CreateSeuratObject(counts = F52G7, min.features = 100)

F52G21 <- Read10X(data.dir = "/work/dimitri/single/2052_matrix/filtered_2052_gex21") 
F52G21_obj <- CreateSeuratObject(counts = F52G21, min.features = 100)

F52G28 <- Read10X(data.dir = "/work/dimitri/single/2052_matrix/filtered_2052_gex28") 
F52G28_obj <- CreateSeuratObject(counts = F52G28, min.features = 100)

# 2053 samples ----
F53G0 <- Read10X(data.dir = "/work/dimitri/single/2053_matrix/filtered_2053_gex0") 
F53G0_obj <- CreateSeuratObject(counts = F53G0, min.features = 100)

F53G7 <- Read10X(data.dir = "/work/dimitri/single/2053_matrix/filtered_2053_gex7") 
F53G7_obj <- CreateSeuratObject(counts = F53G7, min.features = 100)

F53G21 <- Read10X(data.dir = "/work/dimitri/single/2053_matrix/filtered_2053_gex21") 
F53G21_obj <- CreateSeuratObject(counts = F53G21, min.features = 100)

F53G28 <- Read10X(data.dir = "/work/dimitri/single/2053_matrix/filtered_2053_gex28") 
F53G28_obj <- CreateSeuratObject(counts = F53G28, min.features = 100)

# 2055 samples ----
F55G0 <- Read10X(data.dir = "/work/dimitri/single/2055_matrix/filtered_2055_gex0") 
F55G0_obj <- CreateSeuratObject(counts = F55G0, min.features = 100)

F55G7 <- Read10X(data.dir = "/work/dimitri/single/2055_matrix/filtered_2055_gex7") 
F55G7_obj <- CreateSeuratObject(counts = F55G7, min.features = 100)

F55G21 <- Read10X(data.dir = "/work/dimitri/single/2055_matrix/filtered_2055_gex21") 
F55G21_obj <- CreateSeuratObject(counts = F55G21, min.features = 100)

F55G28 <- Read10X(data.dir = "/work/dimitri/single/2055_matrix/filtered_2055_gex28") 
F55G28_obj <- CreateSeuratObject(counts = F55G28, min.features = 100)

# P1 sample ----
FP1G0 <- Read10X(data.dir = "/work/dimitri/single/P1_matrix/filtered_p1_gex0") 
FP1G0_obj <- CreateSeuratObject(counts = FP1G0, min.features = 100)

FP1G7 <- Read10X(data.dir = "/work/dimitri/single/P1_matrix/filtered_p1_gex7") 
FP1G7_obj <- CreateSeuratObject(counts = FP1G7, min.features = 100)

FP1G21 <- Read10X(data.dir = "/work/dimitri/single/P1_matrix/filtered_p1_gex21") 
FP1G21_obj <- CreateSeuratObject(counts = FP1G21, min.features = 100)

FP1G28 <- Read10X(data.dir = "/work/dimitri/single/P1_matrix/filtered_p1_gex28") 
FP1G28_obj <- CreateSeuratObject(counts = FP1G28, min.features = 100)

# P2 sample ----
FP2G0 <- Read10X(data.dir = "/work/dimitri/single/P2_matrix/filtered_p2_gex0") 
FP2G0_obj <- CreateSeuratObject(counts = FP2G0, min.features = 100)

FP2G7 <- Read10X(data.dir = "/work/dimitri/single/P2_matrix/filtered_p2_gex7") 
FP2G7_obj <- CreateSeuratObject(counts = FP2G7, min.features = 100)

FP2G21 <- Read10X(data.dir = "/work/dimitri/single/P2_matrix/filtered_p2_gex21") 
FP2G21_obj <- CreateSeuratObject(counts = FP2G21, min.features = 100)

FP2G28 <- Read10X(data.dir = "/work/dimitri/single/P2_matrix/filtered_p2_gex28") 
FP2G28_obj <- CreateSeuratObject(counts = FP2G28, min.features = 100)

# P3 sample ----
FP3G0 <- Read10X(data.dir = "/work/dimitri/single/P3_matrix/filtered_p3_gex0") 
FP3G0_obj <- CreateSeuratObject(counts = FP3G0, min.features = 100)

FP3G7 <- Read10X(data.dir = "/work/dimitri/single/P3_matrix/filtered_p3_gex7") 
FP3G7_obj <- CreateSeuratObject(counts = FP3G7, min.features = 100)

FP3G21 <- Read10X(data.dir = "/work/dimitri/single/P3_matrix/filtered_p3_gex21") 
FP3G21_obj <- CreateSeuratObject(counts = FP3G21, min.features = 100)

FP3G28 <- Read10X(data.dir = "/work/dimitri/single/P3_matrix/filtered_p3_gex28") 
FP3G28_obj <- CreateSeuratObject(counts = FP3G28, min.features = 100)


# P4 sample ----
FP4G0 <- Read10X(data.dir = "/work/dimitri/single/P4_matrix/filtered_p4_gex0") 
FP4G0_obj <- CreateSeuratObject(counts = FP4G0, min.features = 100)

FP4G7 <- Read10X(data.dir = "/work/dimitri/single/P4_matrix/filtered_p4_gex7") 
FP4G7_obj <- CreateSeuratObject(counts = FP4G7, min.features = 100)

FP4G21 <- Read10X(data.dir = "/work/dimitri/single/P4_matrix/filtered_p4_gex21") 
FP4G21_obj <- CreateSeuratObject(counts = FP4G21, min.features = 100)

FP4G28 <- Read10X(data.dir = "/work/dimitri/single/P4_matrix/filtered_p4_gex28") 
FP4G28_obj <- CreateSeuratObject(counts = FP4G28, min.features = 100)

# P5 sample ----
FP5G0 <- Read10X(data.dir = "/work/dimitri/single/P5_matrix/filtered_p5_gex0") 
FP5G0_obj <- CreateSeuratObject(counts = FP5G0, min.features = 100)

FP5G7 <- Read10X(data.dir = "/work/dimitri/single/P5_matrix/filtered_p5_gex7") 
FP5G7_obj <- CreateSeuratObject(counts = FP5G7, min.features = 100)

FP5G21 <- Read10X(data.dir = "/work/dimitri/single/P5_matrix/filtered_p5_gex21") 
FP5G21_obj <- CreateSeuratObject(counts = FP5G21, min.features = 100)

FP5G28 <- Read10X(data.dir = "/work/dimitri/single/P5_matrix/filtered_p5_gex28") 
FP5G28_obj <- CreateSeuratObject(counts = FP5G28, min.features = 100)

# P6 sample ----
FP6G0 <- Read10X(data.dir = "/work/dimitri/single/P6_matrix/filtered_p6_gex0") 
FP6G0_obj <- CreateSeuratObject(counts = FP6G0, min.features = 100)

FP6G7 <- Read10X(data.dir = "/work/dimitri/single/P6_matrix/filtered_p6_gex7") 
FP6G7_obj <- CreateSeuratObject(counts = FP6G7, min.features = 100)

FP6G21 <- Read10X(data.dir = "/work/dimitri/single/P6_matrix/filtered_p6_gex21") 
FP6G21_obj <- CreateSeuratObject(counts = FP6G21, min.features = 100)

# P7 sample ----
FP7G0 <- Read10X(data.dir = "/work/dimitri/single/P7_matrix/filtered_p7_gex0") 
FP7G0_obj <- CreateSeuratObject(counts = FP7G0, min.features = 100)

FP7G7 <- Read10X(data.dir = "/work/dimitri/single/P7_matrix/filtered_p7_gex7") 
FP7G7_obj <- CreateSeuratObject(counts = FP7G7, min.features = 100)

FP7G21 <- Read10X(data.dir = "/work/dimitri/single/P7_matrix/filtered_p7_gex21") 
FP7G21_obj <- CreateSeuratObject(counts = FP7G21, min.features = 100)

FP7G28 <- Read10X(data.dir = "/work/dimitri/single/P7_matrix/filtered_p7_gex28") 
FP7G28_obj <- CreateSeuratObject(counts = FP7G28, min.features = 100)

# P8 samples ----
FP8G0 <- Read10X(data.dir = "/work/dimitri/single/P8_matrix/filtered_p8_gex0") 
FP8G0_obj <- CreateSeuratObject(counts = FP8G0, min.features = 100)

FP8G7 <- Read10X(data.dir = "/work/dimitri/single/P8_matrix/filtered_p8_gex7") 
FP8G7_obj <- CreateSeuratObject(counts = FP8G7, min.features = 100)

FP8G21 <- Read10X(data.dir = "/work/dimitri/single/P8_matrix/filtered_p8_gex21") 
FP8G21_obj <- CreateSeuratObject(counts = FP8G21, min.features = 100)

FP8G28 <- Read10X(data.dir = "/work/dimitri/single/P8_matrix/filtered_p8_gex28") 
FP8G28_obj <- CreateSeuratObject(counts = FP8G28, min.features = 100)

# P9 samples ----
FP9G0 <- Read10X(data.dir = "/work/dimitri/single/P9_matrix/filtered_p9_gex0") 
FP9G0_obj <- CreateSeuratObject(counts = FP9G0, min.features = 100)

FP9G21 <- Read10X(data.dir = "/work/dimitri/single/P9_matrix/filtered_p9_gex21") 
FP9G21_obj <- CreateSeuratObject(counts = FP9G21, min.features = 100)

FP9G28 <- Read10X(data.dir = "/work/dimitri/single/P9_matrix/filtered_p9_gex28") 
FP9G28_obj <- CreateSeuratObject(counts = FP9G28, min.features = 100)

# Merge the individuals Seurat object into one ----
merged_seurat <- merge(x = F47G0_obj, 
                       y = c(F47G7_obj, F47G21_obj, F47G28_obj, F49G0_obj, F49G7_obj, F49G21_obj, F49G28_obj, F51G0_obj, F51G7_obj, F51G21_obj, F51G28_obj, F52G0_obj, F52G7_obj, F52G21_obj, F52G28_obj, F53G0_obj, F53G7_obj, F53G21_obj, F53G28_obj, F55G0_obj, F55G7_obj, F55G21_obj, F55G28_obj, FP1G0_obj, FP1G7_obj, FP1G21_obj, FP1G28_obj, FP2G0_obj, FP2G7_obj, FP2G21_obj, FP2G28_obj, FP3G0_obj, FP3G7_obj, FP3G21_obj, FP3G28_obj, FP4G0_obj, FP4G7_obj, FP4G21_obj, FP4G28_obj, FP5G0_obj, FP5G7_obj, FP5G21_obj, FP5G28_obj, FP6G0_obj, FP6G7_obj, FP6G21_obj, FP7G0_obj, FP7G7_obj, FP7G21_obj, FP7G28_obj, FP8G0_obj, FP8G7_obj, FP8G21_obj, FP8G28_obj, FP9G0_obj, FP9G21_obj, FP9G28_obj),
                       add.cell.id = c("2047_GEX0","2047_GEX7", "2047_GEX21", "2047_GEX28", "2049_GEX0", "2049_GEX7", "2049_GEX21", "2049_GEX28", "2051_GEX0", "2051_GEX7", "2051_GEX21", "2051_GEX28", "2052_GEX0", "2052_GEX7", "2052_GEX21", "2052_GEX28", "2053_GEX0", "2053_GEX7", "2053_GEX21", "2053_GEX28", "2055_GEX0", "2055_GEX7", "2055_GEX21", "2055_GEX28", "P1_GEX0", "P1_GEX7", "P1_GEX21", "P1_GEX28", "P2_GEX0", "P2_GEX7", "P2_GEX21", "P2_GEX28", "P3_GEX0", "P3_GEX7", "P3_GEX21", "P3_GEX28", "P4_GEX0", "P4_GEX7", "P4_GEX21", "P4_GEX28", "P5_GEX0", "P5_GEX7", "P5_GEX21", "P5_GEX28", "P6_GEX0", "P6_GEX7", "P6_GEX21", "P7_GEX0", "P7_GEX7", "P7_GEX21", "P7_GEX28", "P8_GEX0", "P8_GEX7", "P8_GEX21", "P8_GEX28", "P9_GEX0", "P9_GEX21", "P9_GEX28"))


save(merged_seurat, file="merged_seurat.RData")

