rm *CTA*
rm *png
rm *svg
rm antigen_counts.tsv

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
  gencode.v37.annotation.minimal.gtf | # grep the transcript identifiers from the GTF file
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
  gencode.v37.annotation.minimal.gtf | # grep the transcript identifiers from the GTF file
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
gencode.v37.annotation.minimal.gtf | # grep the transcript identifiers from the GTF file
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
  COORDS=`grep \\"${line}\\" gencode.v37.annotation.minimal.gtf | 
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
    COORDS=`grep \\"${line}\\" gencode.v37.annotation.minimal.gtf | 
    grep "	gene	" | # Only getting gene entries
    cut -f 1,4,5 | # Getting gene coordinates
    sed 's/chr/hs/g'`;  # Convering chr prefix to hs prefix for Circos
    echo "${COORDS}	0" >> ${i}.CTA_coords.tsv; 
  done < ${i}.CTAs ; 
done

echo "Creating Circos plot..."
circos --config nutm1_circos.config
mv circos.svg nutm1_cta_6tracks.svg
mv circos.png nutm1_cta_6tracks.png



echo "Working on barplot..."

echo "Extracting PRAME and NUTM1 transcript identifiers..."

grep '"PRAME"' *gtf |  # Extract PRAME entries
grep "	transcript	" | # Filter for only transcripts
cut -f 9 |  # Get annotation column
cut -f 2 -d ';' |  # Get transcript identifier column
sed 's/ transcript_id "//g' | # Strip off unneeded text
sed 's/"//g' \
> PRAME_tx_ids


grep '"NUTM1"' *gtf |  # Extract NUTM1 entries
grep "	transcript	" |  # filter for only transcripts
cut -f 9 |  # Get annotation column
cut -f 2 -d ';' |  # Get transcript identifier column
sed 's/ transcript_id "//g' | # Strip off unneeded text
sed 's/"//g' \
> NUTM1_tx_ids

cat NUTM1_tx_ids > NUTM1_PRAME_tx_ids
cat PRAME_tx_ids >> NUTM1_PRAME_tx_ids

echo "Generating antigen counts for barplot..."

echo "Line	Antigen_Source	pMHC_Count" > antigen_counts.tsv

# FOR ALL BUT PDX, ALL BUT CTAS
for samp in 797 1015 14169 JCM1 PER403; do  # For each line (except for PDX)
  for i in SPLICE FUSION ERV VIRUS; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep "	${samp}-" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    grep -v '1 ' |  # Remove rows that occur only once (to ensure in at leats 2 reps)
                    wc -l`;  # Get count
    echo ${samp}	${i}	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR PDX, ALL BUT CTAS
for samp in PDX; do  # For each line (except for PDX)
  for i in SPLICE FUSION ERV VIRUS; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep "	${samp}	" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    wc -l`;  # Get count
    echo ${samp}	${i}	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, CTAs except PRAME/NUTM1
for samp in 797 1015 14169 JCM1 PER403; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -v -f NUTM1_PRAME_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}-" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    grep -v '1 ' |  # Remove rows that occur only once (to ensure in at leats 2 reps)
                    wc -l`;  # Get count
    echo ${samp}	${i}-other	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, CTAs except PRAME/NUTM1
for samp in PDX; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -v -f NUTM1_PRAME_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}	" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    wc -l`;  # Get count
    echo ${samp}	${i}-other	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, PRAME
for samp in 797 1015 14169 JCM1 PER403; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -f PRAME_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}-" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    grep -v '1 ' |  # Remove rows that occur only once (to ensure in at leats 2 reps)
                    wc -l`;  # Get count
    echo ${samp}	${i}-PRAME	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, PRAME
for samp in PDX; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -f PRAME_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}	" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    wc -l`;  # Get count
    echo ${samp}	${i}-PRAME	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, NUTM1
for samp in 797 1015 14169 JCM1 PER403; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -f NUTM1_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}-" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    grep -v '1 ' |  # Remove rows that occur only once (to ensure in at leats 2 reps)
                    wc -l`;  # Get count
    echo ${samp}	${i}-NUTM1	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

# FOR ALL BUT PDX, PRAME
for samp in PDX; do  # For each line (except for PDX)
  for i in "CTA/SELF"; do  # For each antigen source
    ANTIGEN_COUNT=`csvcut -d '	' -c peptide,antigen_source,allele,patient_identifier,netmhcpan_4.1b-aff_nm,transcript_id \
                                     nutm1.agg.lens_report.tsv | # Extract relevant columns
                    sed 's/,/	/g' |  # Turn commas to tabs
                    grep -f NUTM1_tx_ids | # Remove CTAs from NUTM1 and PRAME
                    grep "	${samp}	" |  # Grep for samp identifier
                    awk '$5 < 500' | # Filter for binding affinities <500 nM
                    grep "	${i}	" |  # Filter for antigen source of interest
                    cut -f 1,2,3 |  # Extract first three columns
                    sort |  # Sort rows
                    uniq -c |  # Get unique set with count of occurrences
                    wc -l`;  # Get count
    echo ${samp}	${i}-NUTM1	${ANTIGEN_COUNT} >> antigen_counts.tsv; 
  done; 
done

sed -i 's/ /	/g' antigen_counts.tsv

sed -i 's/SPLICE/Splice Variant/g' antigen_counts.tsv
sed -i 's/FUSION/Gene Fusion/g' antigen_counts.tsv
sed -i 's/VIRUS/Virus/g' antigen_counts.tsv
sed -i 's/CTA\/SELF-PRAME/PRAME (CTA)/g' antigen_counts.tsv
sed -i 's/CTA\/SELF-NUTM1/NUTM1 (CTA)/g' antigen_counts.tsv
sed -i 's/CTA\/SELF-other/Other CTA/g' antigen_counts.tsv

sed -i 's/1015/10-15/g' antigen_counts.tsv
sed -i 's/797/TC-797/g' antigen_counts.tsv


/share/apps/clusterapps/R-4.4.1/bin/Rscript make_antigen_barplot.R

echo "Making peptide heatmaps for NUTM1..."

echo "Extracting pMHCs from 797..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv | 
sed 's/,/	/g' | 
awk '$4 ~ /797/' | 
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' |
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/797	/g' > nutm1_matrix_inputs.tsv

echo "Extracting pMHCs from 1015..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /1015/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/1015	/g' >> nutm1_matrix_inputs.tsv

echo "Extracting pMHCs from 14169..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /14169/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/14169	/g' >> nutm1_matrix_inputs.tsv

echo "Extracting pMHCs from PDX..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /PDX/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/PDX	/g' >> nutm1_matrix_inputs.tsv

echo "Extracting pMHCs from JCM1..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /JCM1/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/JCM1	/g' >> nutm1_matrix_inputs.tsv

echo "Extracting pMHCs from PER403..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /PER403/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f NUTM1_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/PER403	/g' >> nutm1_matrix_inputs.tsv

echo "Making peptide heatmaps for PRAME..."

echo "Extracting pMHCs from 797..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv | 
sed 's/,/	/g' | 
awk '$4 ~ /797/' | 
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' |
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/797	/g' > prame_matrix_inputs.tsv

echo "Extracting pMHCs from 1015..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /1015/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/1015	/g' >> prame_matrix_inputs.tsv

echo "Extracting pMHCs from 14169..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /14169/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/14169	/g' >> prame_matrix_inputs.tsv

echo "Extracting pMHCs from PDX..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /PDX/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/PDX	/g' >> prame_matrix_inputs.tsv

echo "Extracting pMHCs from JCM1..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /JCM1/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/JCM1	/g' >> prame_matrix_inputs.tsv

echo "Extracting pMHCs from PER403..."
csvcut -d '	' -c pos,peptide,antigen_source,patient_identifier,transcript_id,netmhcpan_4.1b-aff_nm nutm1.agg.lens_report.tsv |
sed 's/,/	/g' |
awk '$4 ~ /PER403/' |
awk '$6 < 500' | 
grep "CTA/SELF" | 
grep -f PRAME_tx_ids | 
cut -f 1,2,3,5 | 
sort | 
uniq -c  | 
grep -v '1 ' | 
rev | 
cut -f 1 -d ' ' | 
rev | 
sort -n -k1,1 | 
sed 's/^/PER403	/g' >> prame_matrix_inputs.tsv

python3 create_cta_matrix.py -i nutm1_matrix_inputs.tsv -o nutm1_matrix.tsv -l 1132
python3 create_cta_matrix.py -i prame_matrix_inputs.tsv -o prame_matrix.tsv -l 509

sed -i 's/^797/TC-797/g' nutm1_matrix.tsv
sed -i 's/^797/TC-797/g' prame_matrix.tsv

sed -i 's/^1015/10-15/g' nutm1_matrix.tsv
sed -i 's/^1015/10-15/g' prame_matrix.tsv

/share/apps/clusterapps/R-4.4.1/bin/Rscript make_heatmaps.R
