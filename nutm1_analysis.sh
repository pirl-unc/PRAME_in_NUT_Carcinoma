echo "Creating list of CTAs per cell line."

# Get the list of CTAs expressed in each cell line. pMHCs must have <500 nM
# affinity and be observed in at least two replicates for the sample.
for i in 1015 14169 797 JCM1; do 
  echo ${i}
  grep -f <(grep "	${i}\-" nutm1.agg.lens_report.tsv |
            grep "	CTA/SELF	" | # Only CTAs
            awk '$9 < 500' | # Only <500 nM CTAs
            cut -f 30,36 | # Extract transcript identifier and patient identifier
            sort | # sort and unique so that each transcript is only represented once per patient
            uniq | # to prevent double counting due to multiple alleles
            cut -f 1 | # cut transcript identifier
            sort | # sort transcript identifiers
            uniq -c | # Create unique set of transcript identifiers and count occurrences
            grep -v '1 ' | # Remove transcript identifers only observed in a single transcript
            rev | # rev the string
            cut -f 1 -d ' ' | # So that we can easily extract the transcript identifier and rev it
            rev) \
  gencode.v37.annotation.with.hervs.gtf | # grep the transcript identifiers from the GTF file
  cut -f 9 | # Get the annotation column
  cut -f 4 -d ';' | # Get the gene symbol column
  rev | # rev so we can easily...
  cut -f 1 -d ' ' | # cut the gene symbol
  rev | # and rev again so the gene symbol is not backwards
  sed 's/"//g' | # Replace quotes with nothing
  sort | # sort list of gene symbols and list only unique gene symbols
  uniq \
  > ${i}.CTAs; # Output the gene symbol list to <LINE>.CTAs
done

# Get the list of CTAs expressed in PDX. pMHCs must have <500 nM affinity.
# pMHCs cannot be observed in multiple reps since PDX only has a single rep.
# The code has been modified to reflect that.
for i in PDX; do 
  echo ${i}
  grep -f <(grep "	${i}	" nutm1.agg.lens_report.tsv |
            grep "	CTA/SELF	" | # Only CTAs
            awk '$9 < 500' | # Only <500 nM CTAs
            cut -f 30,36 | # Extract transcript identifier and patient identifier
            sort | # sort and unique so that each transcript is only represented once per patient
            uniq | # to prevent double counting due to multiple alleles
            cut -f 1 | # cut transcript identifier
            sort | # sort transcript identifiers
            uniq -c | # Create unique set of transcript identifiers and count occurrences
            rev | # rev the string
            cut -f 1 -d ' ' | # So that we can easily extract the transcript identifier and rev it
            rev) \
  gencode.v37.annotation.with.hervs.gtf | # grep the transcript identifiers from the GTF file
  cut -f 9 | # Get the annotation column
  cut -f 4 -d ';' | # Get the gene symbol column
  rev | # rev so we can easily...
  cut -f 1 -d ' ' | # cut the gene symbol
  rev | # and rev again so the gene symbol is not backwards
  sed 's/"//g' | # Replace quotes with nothing
  sort | # sort list of gene symbols and list only unique gene symbols
  uniq \
  > ${i}.CTAs; # Output the gene symbol list to <LINE>.CTAs
done

# Get the list of CTAs expressed in Per403/PER403. The reps for this sample use
# different capitalization of the sample name, so a unique method is required
echo "PER403/Per403"
grep -f <(grep "R403-" nutm1.agg.lens_report.tsv |
          grep "	CTA/SELF	" | # Only CTAs
          awk '$9 < 500' | # Only <500 nM CTAs
          cut -f 30,36 | # Extract transcript identifier and patient identifier
          sort | # sort and unique so that each transcript is only represented once per patient
          uniq | # to prevent double counting due to multiple alleles
          cut -f 1 | # cut transcript identifier
          sort | # sort transcript identifiers
          uniq -c | # Create unique set of transcript identifiers and count occurrences
          grep -v '1 ' | # Remove transcript identifers only observed in a single transcript
          rev | # rev the string
          cut -f 1 -d ' ' | # So that we can easily extract the transcript identifier and rev it
          rev) \
gencode.v37.annotation.with.hervs.gtf | # grep the transcript identifiers from the GTF file
cut -f 9 | # Get the annotation column
cut -f 4 -d ';' | # Get the gene symbol column
rev | # rev so we can easily...
cut -f 1 -d ' ' | # cut the gene symbol
rev | # and rev again so the gene symbol is not backwards
sed 's/"//g' | # Replace quotes with nothing
sort | # sort list of gene symbols and list only unique gene symbols
uniq \
> PER403.CTAs; # Output the gene symbol list to <LINE>.CTAs

echo "Creating union of CTAs observed across cell lines..."
cat *.CTAs | sort | uniq > all.CTAs

echo "Creating CTA label track..."
while read line; 
  do echo ${line}; 
  COORDS=`grep \\"${line}\\" gencode.v37.annotation.with.hervs.gtf | 
  grep "	gene	" | # Only getting gene entries
  cut -f 1,4,5 | # Getting gene coordinates
  sed 's/chr/hs/g'`;  # Convering chr prefix to hs prefix for Circos
  echo "${COORDS}	${line}" >> all_CTA_coords.tsv; 
done < all.CTAs

echo "Creating Circos tracks for each cell line..."
for i in 1015 14169 797 PDX JCM1 PER403; do
  echo ${i}
  while read line; 
    do echo ${line}; 
    COORDS=`grep \\"${line}\\" gencode.v37.annotation.with.hervs.gtf | 
    grep "	gene	" | # Only getting gene entries
    cut -f 1,4,5 | # Getting gene coordinates
    sed 's/chr/hs/g'`;  # Convering chr prefix to hs prefix for Circos
    echo "${COORDS}	0" >> ${i}.CTA_coords.tsv; 
  done < ${i}.CTAs ; 
done
