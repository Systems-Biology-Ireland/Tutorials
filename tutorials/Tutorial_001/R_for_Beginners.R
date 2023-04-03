# Tutorial for R
# Luis Iglesias-Martinez

# This script is meant to introduce concepts such as variables, classes, functions and variables. 
# It is an introduction to R.

# Lines that start with # are called comments and they don't run in the console. 

# We will start with creating a nummeric variable called A with value of 0
# You can copy each line into the console 

A = 0

# Or click control enter.

class(A)

# This will print the type of the variable. 

# Check in your environment and see that now you have a variable called A with value of 0.

# Now we will showcase some basic arithmetic operations with R, like the addition:

A + 1

# The result is printed underneath. It is 1. But see that we haven't stored this value anywhere.
# Now, let's create a variable to store it:

B = A  + 1

# You can also change the values of a variable

A = B + A 

# Check the environment to the left. Now you have variable B with value of 1.

# There's more operations you can do 

A - B

# And multiplications or divisions

B*A

A/B

# We can also use () to create more complex arithmetic operations.

C = B*(A + 1)

# And ^ can be used for power

D = C^2

# And log() for logarithms

E = log(D)

# exp for exponential.

E = exp(E)

# You can also create vectors.

F = c(0,1,0)

# They are objects with several entries. The entries can be accessed with []

F[1]

F[3]

# if you try to access an entry that doesn't exist, you will get an NA.

F[4]

# Now let's create a matrix. We can use the function matrix(value, dimension1, dimension2)

G = matrix(0,3,3)

# A matrix has rows and columns that you can access using matrix[row,column] as follows:

G[1,1]

# that is the value in the  first row and the first column.

# You can change values in the matrix too.

G[1,1] = 1

# And you can access complete row or complete column.

G[1,] # accessing first row.

G[,2] # accessing second column.

# and change their values. 

G[2,] = 2

# You can check if values are the equal, lower or greater using logic operations.

G[1,1] == 1

G[1,1] <  0

G[1,1] >= 0

# This are called boolean data. You can also store them in a variable.

H = G[1,1] == 1

# You can also transform it to 1 and 0. (1 for TRUE) By multiplying by 1.

1*H

# Now let's create a string or a charater.

H = 'a'

# This means that the value of H is now "a". See how you can't combines data of numeric and character type.

H + A

# This outputs an error as you saw.

# You can combine strings with strings though.

I = "b"

J = paste(H,I, sep = "")

# This is useful to find files.

# We will now download the gene expression data of uveal melanoma patients from TCGA.

# You can download it from the Broad's Institute Firehose website.

# You can also download it your machine using R.

url_uvm = 'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/UVM/20160128/gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz'


download.file(url_uvm, 'uvm_rnaseq.tar.gz')


# If  this doesn't work for you, you can downaload into downloads and then put the path in untar.

dir.create("uvm_tcga")

# untar 
 
untar('uvm_rnaseq.tar.gz', exdir = "uvm_tcga")

# Now we will read the file we need.
# The getwd() will give you the current working directory. We will paste it to the new folder we've created. 

path_to_file = paste(getwd(),'/uvm_tcga/gdac.broadinstitute.org_UVM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep = "")

# read into workspace. 

uvm = read.delim(path_to_file)

# Now we can inspect it. 

colnames(uvm)

uvm[1,] 

unique(as.vector(uvm[1,]))

# Next, we can focus only on the columns with raw gene counts. 

# First let's store the gene names.

gene_ids = uvm[,1]


raw_counts = uvm[1,] == 'raw_count'

# we can select the columns with raw counts. 


uvm_rc = uvm[,raw_counts]

# We will delete the first row as we don't need it in uvm_rc

uvm_rc = uvm_rc[-1,]

# we will do the same for the gene_ids

gene_ids = gene_ids[-1]

# give the dataframe genes rownames

rownames(uvm_rc) = as.vector(gene_ids)


# The entries of the dataframe are strings. 

uvm_rc[1,1]

# we need to make them numbers.

uvm_rc = mapply(as.numeric, uvm_rc)

# We can make a simple plot of the first row. 

plot(uvm_rc[1000,])

# Now we can also store the means of each sample and the variance.

uvm_means = matrix(0,dim(uvm_rc)[2],1)
uvm_variances = matrix(0,dim(uvm_rc)[2],1)

# There's other ways of doing this but I want to introduce another topic

# For loops.
for (i in 1:dim(uvm_rc)[2]){
  uvm_means[i] = mean(uvm_rc[,i])
  uvm_variances[i] = var(uvm_rc[,i])
}

# Let's plot a histogram.

hist(uvm_variances)

# this shows how the variance is spread in each patient.

# Let's download Bioconductor.
install.packages('BiocManager') 

BiocManager::install('DESeq2')

# Load pacakge to use functions. 

library(DESeq2)

# create a dds object.
# we're rounding the counts as they need to be integers. 

coldata = matrix(0,80,1)

colnames(coldata) = 'cols'# just a random name to use itin DESeq2.

# We will make a random design as for this, we are not looking to do differential expression analysis. 
# However, you can change this and make a real design and coldata when doing differential expression analysis.

dds  = DESeqDataSetFromMatrix(countData = round(uvm_rc),colData = coldata, design = ~1)

# now we can now normalize.

dds =  DESeq(dds)

# Get the variance stabilization transformed data.

vst =  vst(dds, blind=FALSE)

# We can get the normalized values through

vst_a = assay(vst)

# Calculate the new variances.



uvm_vst_means = matrix(0,dim(uvm_rc)[2],1)
uvm_vst_variances = matrix(0,dim(uvm_rc)[2],1)

# There's other ways of doing this but I want to introduce another topic

# For loops.
for (i in 1:dim(uvm_rc)[2]){
  uvm_vst_means[i] = mean(vst_a[,i])
  uvm_vst_variances[i] = var(vst_a[,i])
}

# check histogram

hist(vst)
