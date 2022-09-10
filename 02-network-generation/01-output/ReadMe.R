freshwater = read.table('freshwater.pa.tsv', sep='\t', row.names=1, header=1)
hypothermal.nonmarine = read.table('hypothermal-nonmarine.pa.tsv', sep='\t', row.names=1, header=1)
marine.water = read.table('marine-water.pa.tsv', sep='\t', row.names=1, header=1)
saline.hypersaline = read.table('saline-hypersaline.pa.tsv', sep='\t', row.names=1, header=1)
thermal = read.table('thermal.pa.tsv', sep='\t', row.names=1, header=1)
host.associated = read.table('host-associated.pa.tsv', sep='\t', row.names=1, header=1)
marine.sediment = read.table('marine-sediment.pa.tsv', sep='\t', row.names=1, header=1)
oil = read.table('oil.pa.tsv', sep='\t', row.names=1, header=1)
soils = read.table('soils.pa.tsv', sep='\t', row.names=1, header=1)
water.treatment = read.table('water-treatment.pa.tsv', sep='\t', row.names=1, header=1)

dim(freshwater)
dim(hypothermal.nonmarine)
dim(marine.water)
dim(saline.hypersaline)
dim(thermal)
dim(host.associated)
dim(marine.sediment)
dim(oil)
dim(soils)
dim(water.treatment)

length(unique(c(rownames(freshwater), rownames(hypothermal.nonmarine), rownames(marine.water), rownames(saline.hypersaline), 
                rownames(thermal), rownames(host.associated), rownames(marine.sediment), rownames(oil), rownames(soils), rownames(water.treatment))))


