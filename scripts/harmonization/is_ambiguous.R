# Function to identify ambiguous SNPs
is_ambiguous <- function(allele1, allele2) {
  (allele1 == "A" & allele2 == "T") | (allele1 == "T" & allele2 == "A") |
    (allele1 == "G" & allele2 == "C") | (allele1 == "C" & allele2 == "G")
}

