# This file contains notes of command-line work for the Peanut Core Collection Genotyping Project (PCCGP)
# of Steven Cannon, Ethy Cannon et al., 2019-2020, for the manuscript
# "Genotypic characterization of the U.S. peanut core collection" (Otyama, Kulkarni et al., in preparation early 2020).

# The notes are mostly chronological. Although they faithfully record the command-line data explorations,
# they are not particularly well documented. Reading the notes will require some fluency with 
# Unix and computational biology. If you have questions about particular steps, please
# contact Steven Cannon: scannon@iastate.edu or steven.cannon@usda.gov


# Started 2019-10-05

# Note regarding identities of four samples with suspicious PI labels, 2019-10-28:
# After sleuthing recorded in case_of_the_missing_PIs.xlsx, identified 3 accessions that were mis-transcribed,
# with positive identification of the accessions that these should have been,
# and a fourth that was mis-transcribed but with uncertain identification as to the correct accession.
# Swapped accessions PI345401, PI442697, PI473562 and removed PI494595 from all files. Re-calculated analyses as needed.
  # The incorrectly identified accessions:
  #   PI 345401 - Triticum
  #   PI 442697 - Solanum
  #   PI 473562 - Rice
  #   PI 494595 - Sunflower
  # 
  #              should be:
  #   PI345401 - PI343401_3  red
  #   PI442697 - PI442597_2  blue
  #   PI473562 - PI493562_3  green
  #   PI494595 - PI493595??  - but REMOVE from analyses because this one is uncertain (only 91% identity)
  # 

#####
# Find alleles corresponding with the Affy SNPs, from the diploid genomes (A. ipaensis and A. duranensis).
#
# Get flanking sequence from SNPs_whole.txt, and use this to search the diploid genomes.
# Search the v2 assemblies, with the reasoning that those alleles should be at least as "correct" as in the v1 assemblies.
# Make two versions: one with the first SNP used in the sequence, and one with the second SNP used in the sequence.

  mkdir 02_genomic_SNPs

  # With both SNPs reported, in [A/G] form:
    cat 00_start/SNPs_whole.txt | awk '$1~/^AX/ {print ">" $1 " " $1123 " " $1124 "\n" $1126}' |
      cat > 02_genomic_SNPs/48k_SNPs_flanking.fna
  # With the first SNP incorporated:
    cat 00_start/SNPs_whole.txt | awk -v OFS="\t" '$1~/^AX/ {print $1, $1123, $1124, $1126}' |
      perl -pe 's/\[(\w)\/(\w)\]/\t$1\t$2\t/' |
      awk -v OFS="\t" '{print $1 "." $5, $2, $3, $4 $5 $7}' > 02_genomic_SNPs/48k_SNPs_flanking_1stSNP.fas0
  # With the second SNP incorporated:
    cat 00_start/SNPs_whole.txt | awk -v OFS="\t" '$1~/^AX/ {print $1, $1123, $1124, $1126}' |
      perl -pe 's/\[(\w)\/(\w)\]/\t$1\t$2\t/' |
      awk -v OFS="\t" '{print $1 "." $6, $2, $3, $4 $6 $7}' > 02_genomic_SNPs/48k_SNPs_flanking_2ndSNP.fas0
  # Combine 1st and 2nd:
    cat 00_start/48k_SNPs_flanking_1stSNP.fas0 02_genomic_SNPs/48k_SNPs_flanking_2ndSNP.fas0 | sort |
      awk '{print ">" $1 " " $2 " " $3 "\n" $4}' > 02_genomic_SNPs/48k_SNPs_combined.fna

  mkdir 01_dip_genomes
  cd 01_dip_genomes
    ln -s ../../araip_dovetail/K30076.gnm2.1GWY/K30076.gnm2.1GWY.genome_main.fna
    ln -s ../../aradu_dovetail/V14167.gnm2.J7QH/V14167.gnm2.J7QH.genome_main.fna
    cd ..

 # Get these genome assemblies from lis-stage
  mkdir 01_tetraploids
  ls 01_tetraploids
    arahy.Fuhuasheng.gnm1.XX5Y.genome_main.fna
    arahy.Shitouqi.gnm1.L4VP.genome_main.fna
    arahy.Tifrunner.gnm2.J5K5.genome_main.fna

  mkdir blastdb blastout

  for path in 01_dip_genomes/*fna; do
    base=`basename $path .genome_main.fna`
    echo $base
    makeblastdb -in $path -dbtype nucl -hash_index -parse_seqids -title $base -out blastdb/$base &
  done

  for path in 01_genomes/*fna; do
    base=`basename $path .genome_main.fna`
    echo $base
    makeblastdb -in $path -dbtype nucl -hash_index -parse_seqids -title $base -out blastdb/$base &
  done


  qry_base=48k_SNPs_combined
  sbj_base=K30076.gnm2.1GWY
  blastn -query 02_genomic_SNPs/$qry_base.fna \
    -db blastdb/$sbj_base -perc_identity 95 -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -num_threads 6 -out blastout/$qry_base.x.$sbj_base.bln &

  qry_base=48k_SNPs_combined
  sbj_base=V14167.gnm2.J7QH
  blastn -query 02_genomic_SNPs/$qry_base.fna \
    -db blastdb/$sbj_base -perc_identity 95 -evalue 1e-10 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -num_threads 6 -out blastout/$qry_base.x.$sbj_base.bln &

  qry_base=48k_SNPs_combined
  for path in 01_genomes/*fna; do
    sbj_base=`basename $path .genome_main.fna`
    blastn -query 02_genomic_SNPs/$qry_base.fna \
      -db blastdb/$sbj_base -perc_identity 95 -evalue 1e-10 \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
      -num_threads 6 -out blastout/$qry_base.x.$sbj_base.bln &
  done

##########
# Reshape and filter the BLAST output of the flanking sequences. 
# These had two query sequences per SNP location: one for each variant.
# Suppress these markers, which are repetitive: AX-176806741, AX-176810276, AX-176807409
# Suppress scaffold/unplaced markers.
# Modify chromosome names to all have A01 / B01 form.

  mkdir 02_SNP_matches

  cat blastout/48k_SNPs_combined.x.V14167.gnm2.J7QH.bln | 
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/; s/gnm2.Aradu/gnm2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/aradu.V14167.gnm2.//' |
      awk '$1!~/176806741|176810276|176807409/' |
      sort -k1,1 -k3,3 | top_line.awk| cut -f1,2,3,13 > 02_SNP_matches/SNPs.dip_dur.A

  cat blastout/48k_SNPs_combined.x.K30076.gnm2.1GWY.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/araip.K30076.gnm2.//' |
      awk '$1!~/176806741|176810276|176807409/' | 
      sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.dip_ipa.B



  cat blastout/48k_SNPs_combined.x.arahy.Fuhuasheng.gnm1.XX5Y.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/arahy.Fuhuasheng.gnm1.//' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^A/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Fuh.A

  cat blastout/48k_SNPs_combined.x.arahy.Fuhuasheng.gnm1.XX5Y.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/arahy.Fuhuasheng.gnm1.//' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^B/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Fuh.B



  cat blastout/48k_SNPs_combined.x.arahy.Shitouqi.gnm1.L4VP.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/B1/B0/; s/B20/B10/' | perl -pe 's/arahy.Shitoqi.gnm1.//' | awk '$3!~/Unplaced/' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^A/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Shi.A

  cat blastout/48k_SNPs_combined.x.arahy.Shitouqi.gnm1.L4VP.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/B1/B0/; s/B20/B10/' | perl -pe 's/arahy.Shitoqi.gnm1.//' | awk '$3!~/Unplaced/' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^B/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Shi.B



  cat blastout/48k_SNPs_combined.x.arahy.Tifrunner.gnm2.J5K5.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/\S+Arahy.10/A10/; s/\S+Arahy.0/A0/; s/\S+Arahy.1/B0/; s/\S+Arahy.20/B10/' | awk '$3!~/scaff/' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^A/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Tif.A

  cat blastout/48k_SNPs_combined.x.arahy.Tifrunner.gnm2.J5K5.bln |
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/' | sort -k1,1 -k4nr,4nr |
    awk -v OFS="\t" '$5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, substr($14,36,1), $10, $11 }' |
      perl -pe 's/\S+Arahy.10/A10/; s/\S+Arahy.0/A0/; s/\S+Arahy.1/B0/; s/\S+Arahy.20/B10/' | awk '$3!~/scaff/' |
      awk '$1!~/176806741|176810276|176807409/' |
      awk '$3~/^B/' | sort -k1,1 -k3,3 | top_line.awk | cut -f1,2,3,13 > 02_SNP_matches/SNPs.Tif.B



# Join the basic results for all genomes. Join onto a backbone of all unique SNP+chr IDs.
  mkdir 02_SNP_matches_joined 03_SNPs_four_genomes

  cat 02_SNP_matches/SNPs* | cut -f1 | sort -u > lis.SNP_chr_unq
  
  for path in 02_SNP_matches/*; do
    file=`basename $path`
    join -a1 lis.SNP_chr_unq $path | sed 's/ /\t/g' |
      awk -v OFS="\t" 'NF==4 {print $0} NF==1 {print $1, "___\t___\t___"}' > 02_SNP_matches_joined/$file 
  done

  # Also add the diploid SNPs from the original work in notes_peanit_SNPs_v##.sh
    join -a1 lis.SNP_chr_unq 02_genomic_SNPs/V14167.gnm2.J7QH.snps |
      awk -v OFS="\t" 'NF==2 {print $1, $2} NF==1 {print $1, "___"}' > 02_SNP_matches_joined/SNPs1.dur_V14167.tsv

    join -a1 lis.SNP_chr_unq 02_genomic_SNPs/K30076.gnm2.1GWY.snps |
      awk -v OFS="\t" 'NF==2 {print $1, $2} NF==1 {print $1, "___"}' > 02_SNP_matches_joined/SNPs1.ipa_K30076.tsv


  paste 02_SNP_matches_joined/* > 03_SNPs_four_genomes/SNPs_four_genomes.tsv
  # Manually (in Excel: SNP_QC_Jan25.xlsx) rearrange columns to get these as the first columns:
  # #marker  dip_dur.A  Fuh.A  Shi.A  Tif.A  dip_ipa.B  Fuh.B  Shi.B  Tif.B



#####
# Extract SNPs from Affy main results file, SNPs_whole.txt

# Some testing, to understand structure of file SNPs_whole.txt
    cat 00_start/SNPs_whole.txt | awk -v FS="\t" '{print $1134 "\t" $1135 "\t" $1136}' | 
      perl -pe 's/\w+\[(\w)\/(\w)\]\w+/$1\t$2/'
       # Flank  Allele_A  Allele_B
       # C  T  T  C
       # G  T  T  G
       # A  G  A  G
       # C  T  T  C
       
    cat 00_start/SNPs_whole.txt | awk -v FS="\t" '{print $1134 "\t" $1135 "\t" $1136}' | 
      perl -pe 's/\w+\[(\w)\/(\w)\]\w+/$1\t$2/' |
      awk '{print "AA:" $1 $1 "\tBB:" $2 $2}'
       # AA:CC  BB:TT
       # AA:GG  BB:TT
       # AA:AA  BB:GG
       # AA:CC  BB:TT

# Convert from AA/BB calls to ATCG calls, from SNPs_whole.txt

    cat 00_start/SNPs_whole.txt | 
      perl -a -F"\t" -ne '
        if ($F[0]=~/^#/){next}
        elsif ($F[0]=~/probeset/){
          $header=$_;
          $header=~s/.CEL_call_code//g;
          $header=~s/\tAffy_SNP_ID.+//;
          print $header;
        }
        elsif ($F[0]=~/^AX/){
          $A=$F[1134]; $B=$F[1135]; 
          for($i=0; $i<=1120; $i++){
            if ($i==0){print "$F[$i]\t"}
            else {
              if    ($F[$i]=~/AA|AB/) {print "$A\t"} 
              elsif ($F[$i]=~/BB/)    {print "$B\t"} 
              elsif ($F[$i]=~/NoCall/){print "-\t"} 
              else                    {print "X\t"}
            }
          }
          print "\n"; 
        } ' | perl -pe 's/^probeset/#probeset/; s/\t+$//' | sort -k1,1 > 00_start/SNPs_alleles_rowsSNPs.tsv 


# Note 2020-01-26: As of this date, I have re-done all steps to this point, in a new, clean
# work directory: peanut_core ... moving the previous one to peanut_coreBAK .
# The reason for the rework is to correct incongruence between the extracted SNPs from the genome assemblies
# and the SNP matrix.
        
# Join the Affy SNPs from the core ... to the SNPs from the diploid assemblies
# 2020-01-26  Use SNPs from four genomes, in SNPs_four_genomes.tsv; also see notes at notes_check_tetraploid_SNPs_v04.sh

# I see in peanut_core_v11.xlsx "notes": 
#   Nov 18: removed PI493562_1, because seems to be artifactual derivative of PI493562_2 and PI493562_3
#   That is a550846-4392289-041321-474_A18; remove row after transposing.
# TO DO: remove PI494595


  # Manual selection of columns in SNP_QC_Jan25.xlsx
  # giving file of SNPs SNPs_four_genomes_to_join.tsv
    sort -o 03_SNPs_four_genomes/SNPs_four_genomes_to_join.tsv 03_SNPs_four_genomes/SNPs_four_genomes_to_join.tsv

    join -a1 00_start/SNPs_alleles_rowsSNPs.tsv 03_SNPs_four_genomes/SNPs_four_genomes_to_join.tsv | 
      perl -pe 's/ /\t/g; s/___/-/g' | # Next: if no SNPs from genomes, report dashes
      awk -v OFS="\t" 'NF==1121 {print $0, "-\t-\t-\t-\t-\t-\t-\t-"} NF==1138 {print $0}' |
      transpose.pl | 
      grep -v "a550846-4392289-041321-474_A18" | sort | awk '$1!~/^PI/' > tmp.SNPs_w_4_genomes

  # hash in the accession IDs, replacing Affy IDs with PI numbers. Move ID to the first column.
    join  tmp.SNPs_w_4_genomes 00_start/ID_peanut_core_v11.hsh | perl -pe 's/ /\t/g' |
      perl -pe 's/^(\S+)\t(.+)\t(\S+)$/$3\t$2/' > 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv

  # Move aa_probeset_id to the first row
    grep aa_probeset_id 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv > tmp.header
    grep -v aa_probeset_id 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv > tmp.rest
    cat tmp.header tmp.rest > 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv

# MAJOR PRODUCT: matrix of SNPs, including from duranensis and ipaensis (in the last column).
#  >>>>> 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv

# Also make a version with SNPs from four genome sequences

# Extract SNPs representative of the clades identified earlier:

# X     aa_probeset_id
# 2     Tif
# 1.1   Fuh
# 1.1   Shi
# 4     dip_dur
# 4     dip_ipa
# 4     PI476025_1
# 4     PI497631_2
# 3     PI468213_s
# 3     PI475971_1
# 3     PI270923_2
# 2     PI210829_2
# 2     PI502120_s
# 1.1   PI268659_1
# 1.1   PI497610_s
# 1.1   PI505975_s
# 1     PI493562_3
# 1     PI341113_2

# Make a tmp file with these 14 names and grep them from SNPs_w_4_genomes.tsv

  grep -f tmp.accns_to_extract 03_SNPs_four_genomes/SNPs_w_4_genomes.tsv |
    transpose.pl > 03_SNPs_four_genomes/SNPs_w_4_genomes_clade_samples_vert.tsv

# Put into spreadsheet chip_and_genome_samples_v02.xlsx
# Merge SNPs from genome samples A and B subgenomes, so A allele if present, otherwise B, otherwise dash
# Filter on genome loci to no-gaps. Copy to text file SNPs_four_genomes_to_join.tsv, then join with main matrix

  # Join SNPs from the genomes onto the main SNP matrix
    join 03_SNPs_four_genomes/SNPs_four_genomes_to_join.tsv 00_start/SNPs_alleles_rowsSNPs.tsv |
      perl -pe 's/ /\t/g' | transpose.pl | awk '$1!~/probeset/ && $1!~/^PI/' | sort -k1,1 > tmp.SNPs_w_4_genomes_merged

  # hash in the accession IDs, replacing Affy IDs with PI numbers. Move ID to the first column.
    join  tmp.SNPs_w_4_genomes_merged 00_start/ID_peanut_core_regions_v11.hsh | perl -pe 's/ /\t/g' |
      grep -v REMOVE |
      perl -pe 's/^(\S+)\t(.+)\t(\S+)$/$3\t$2/' > 03_SNPs_four_genomes/SNPs_w_4_gnm_mrgd.fas

  # Convert to fasta
    cat 03_SNPs_four_genomes/SNPs_w_4_gnm_mrgd.fas |
      perl -pe 's/^(\S+)\t/>$1\n/; s/\t//g' > 04_fas/SNPs_w_4_gnm_mrgd.fas
      # This gives an alignment with 10278 positions and 1123 sequences
 
# Calculate trees
  mkdir 05_trees

  nohup FastTreeMP -nt -quiet 04_fas/SNPs_w_4_gnm_mrgd.fas > 05_trees/SNPs_w_4_gnm_mrgd.tree &

# In archaeopteryx, root tree, and save out in newick format: 05_trees/SNPs_w_4_gnm_mrgd_rt.nh

# Extract tree order
  cat 05_trees/SNPs_w_4_gnm_mrgd_rt2.nh | 
    perl -pe 's/[,;)(]/\n/g' | perl -pe 's/^\d\.\d+:\d.+//; s/:\d\.\d.+//' | 
    sed '/^$/d' > 05_trees/SNPs_w_4_gnm_mrgd_rtr2.tree_order


  # Collapse higly similar sequences
    # Make version without dashes
      cp 04_fas/SNPs_w_4_gnm_mrgd.fas 04_fas/SNPs_w_4_gnm_mrgd_dash.fas
      perl -pi -e 's/-/n/g' 04_fas/SNPs_w_4_gnm_mrgd.fas

    nohup vsearch --cluster_fast 04_fas/SNPs_w_4_gnm_mrgd.fas --id 0.99 --fasta_width 0 \
      --centroids 04_fas/SNPs_w_4_gnm_mrgd_cen99.fas \
      --consout 04_fas/SNPs_w_4_gnm_mrgd_cons99.fas --clusterout_id \
      --uc 04_fas/SNPs_w_4_gnm_mrgd_cons99.uc &

    # Result:
      Clusters: 680 Size min 1, max 144, avg 1.7
      Singletons: 566, 50.4% of seqs, 83.2% of clusters
    
    # clean up centroid deflines
      perl -pi -e 's/>centroid=([^;]+);/>$1 centroid;/' 04_fas/SNPs_w_4_gnm_mrgd_cons99.fas
  
# Calculate tree for reduced alignment  
  nohup FastTreeMP -nt -quiet 04_fas/SNPs_w_4_gnm_mrgd_cen99.fas > 05_trees/SNPs_w_4_gnm_mrgd_cen99.tree &
  nohup FastTreeMP -nt -quiet 04_fas/SNPs_w_4_gnm_mrgd_cons99.fas > 05_trees/SNPs_w_4_gnm_mrgd_cons99.tree &


# Also cluster at 98% identity, for Structure visualization
    nohup vsearch --cluster_fast 04_fas/SNPs_w_4_gnm_mrgd.fas --id 0.98 --fasta_width 0 \
      --centroids 04_fas/SNPs_w_4_gnm_mrgd_cen98.fas \
      --consout 04_fas/SNPs_w_4_gnm_mrgd_cons98.fas --clusterout_id \
      --uc 04_fas/SNPs_w_4_gnm_mrgd_cons98.uc &

    # Result:
      Clusters: 518 Size min 1, max 185, avg 2.2
      Singletons: 387, 34.5% of seqs, 74.7% of clusters

    # clean up centroid deflines
      perl -pi -e 's/>centroid=([^;]+);/>$1 centroid;/' 04_fas/SNPs_w_4_gnm_mrgd_cons98.fas

    # Get list of accessions, for use in Structure
      grep '>' 04_fas/SNPs_w_4_gnm_mrgd_cons98.fas | perl -pe 's/>(\w+)__.+/$1/' > 04_fas/SNPs_w_4_gnm_mrgd_cons98_list_plain
    # Subset at 98% identity, in tree-order:
      grep -f 04_fas/SNPs_w_4_gnm_mrgd_cons98_list_plain 05_trees/SNPs_w_4_gnm_mrgd_rt3.tree_order |
        cat > 05_trees/SNPs_w_4_gnm_mrgd_rt3.tree_order_sub_at_98pct
  

# In Archaeopteryx, root on the clade that contains duranensis_and_ipaensis, then extract tree order:
  cat 05_trees/SNPs_w_4_gnm_mrgd.pxml | grep name |
    perl -pe 's/^ +<name>([^<]+)<.name>/$1/' | sed '/^[01]/d' > 05_trees/SNPs_w_4_gnm_mrgd.tree_order.lis

I AM HERE

#####
# Calculate similarities between all accessions - fasta with 14430 characters

# First, simplify IDs, to simplify later analysis
  hash_into_fasta_id.pl -fas 04_fas/SNPs_alleles_w_dip_simple1.fas -hash 04_fas/hsh.IDs_simplify_dups -out 04_fas/SNPs_alleles_w_dip_simple2.fas
  

  path=04_fas/SNPs_alleles_w_dip_simple2.fas
  base=`basename $path .fas`
  echo $base
  makeblastdb -in $path -dbtype nucl -hash_index -parse_seqids -title $base -out blastdb/$base

  qry_base=SNPs_alleles_w_dip_simple2
  sbj_base=SNPs_alleles_w_dip_simple2
  nohup blastn -query 04_fas/$qry_base.fas \
    -db blastdb/$sbj_base -evalue 1e-10 -outfmt 6 \
    -num_threads 12 -out blastout/$qry_base.x.$sbj_base.bln &



cat blastout/SNPs_alleles_w_dip_simple2.x.SNPs_alleles_w_dip_simple2.bln | awk '$4>14000 {print $3}' |cut -f1 -d'.'| histogram -n -s1 | histplot -d 200
bin abcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRSTabcdefghijKLMNOPQRST
70.00 
71.00 
72.00 
73.00 .
74.00 ...
75.00 ....
76.00 ......
77.00 .............
78.00 ...................
79.00 .........
80.00 ............
81.00 ...........
82.00 .......................
83.00 .............................
84.00 .................................
85.00 ...................
86.00 .....................
87.00 .........................
88.00 .............................
89.00 ....................................
90.00 ..........................................
91.00 .................................
92.00 ......................
93.00 .................
94.00 ..................
95.00 ......................
96.00 ..........................
97.00 .................................................
98.00 ....................................................................
99.00 .........................................................................
100.00  ..


# Cluster at several identity levels, with single-linkage clustering
  mkdir 07_sl_clust

  # for 14430 characters
  cat blastout/SNPs_alleles_w_dip_simple2.x.SNPs_alleles_w_dip_simple2.bln | awk '$3>=99 {print $1 "\t" $2}' > 07_sl_clust/pairs_full.id99.tsv
  blinkPerl_v1.1.pl -in 07_sl_clust/pairs_full.id99.tsv -out 07_sl_clust/clust_full.id99

  cat blastout/SNPs_alleles_w_dip_simple2.x.SNPs_alleles_w_dip_simple2.bln | awk '$3>=98 {print $1 "\t" $2}' > 07_sl_clust/pairs_full.id98.tsv
  blinkPerl_v1.1.pl -in 07_sl_clust/pairs_full.id98.tsv -out 07_sl_clust/clust_full.id98

# And get (non-self) top match. Link the first two fields into one, prior to joining with possible pairs.
# I am doing some extra work to find all matches between PIs-with-replicates. Secondary matches for a given pair are suppressed.
  perl -pe 's/^(\S+)_(\w)\t(\S+)_(\w)/$1\t$2\t$3\t$4/' blastout/SNPs_alleles_w_dip_simple2.x.SNPs_alleles_w_dip_simple2.bln | 
    awk '$1==$3 && $2!=$4' | perl -pe 's/(\S+)\t(\w)\t(\S+)\t(\w)\t/$1_$2__$3_$4\t/' | top_line.awk |
    sort -k1,1 > 07_sl_clust/SNPs_alleles_w_dip_simple2.x.self.all_rep_matches


# Get set of all possible pairings among replicates
  grep '>' 04_fas/SNPs_alleles_w_dip_simple2.fas | sed 's/>//; /_s/d' | sort > tmp.all_rep_IDs1
  grep '>' 04_fas/SNPs_alleles_w_dip_simple2.fas | sed 's/>//; /_s/d' | sort > tmp.all_rep_IDs2

  join -1 99 -2 99 tmp.all_rep_IDs1 tmp.all_rep_IDs2 | perl -pe 's/^ (\S+)_(\w) (\S+)_(\w)$/$1\t$2\t$3\t$4/' | 
    awk '$1==$3 && $2<$4 {print $1 "_" $2 "__" $3 "_" $4}' | sort > all_rep_pairs.tsv
  # looks like this: 
    PI152105_1__PI152105_2
    PI152105_1__PI152105_3
    PI152105_1__PI152105_4
    PI152105_2__PI152105_3
    PI152105_2__PI152105_4
    PI152105_3__PI152105_4

# Join the set of all possible pairings of replicates and the top (non-self) BLAST results
  join -a 1 all_rep_pairs.tsv 07_sl_clust/SNPs_alleles_w_dip_simple2.x.self.all_rep_matches |
    perl -pe 's/__/\t/; s/ /\t/g' > all_rep_pairs_and_all_rep_matches.tsv

# Stats: 
  cat all_rep_pairs_and_all_rep_matches.tsv | awk 'BEGIN{print "matches nomatch"} NF==12 {matches++} NF==2 {nomatch++} END{print matches, nomatch}' 
    matches nomatch
    370 60

##########
# Dec 25
# Check similarity between duranensis_and_ipaensis and other accessions.

# Need to find SNP positions
#   00_start/SNPs_alleles_rowsSNPs.tsv
# Get the positions from the diploid assembly, associated with SNPs ... from BLAST results
  top_line.awk blastout/48k_SNPs_combined.x.* | cut -f1,2,9,10 | awk '{printf "%s\t%s\t%i\n", $1, $2, ($3+$4)/2}' > tmp.48k_SNPs_combined.posns

# Get chromosomal positions from SNPs_whole.txt
  grep -v "^#" 00_start/SNPs_whole.txt | cut -f1,1123,1124 | perl -pe 's/probeset_id/SNP_ID/' | sort > 00_start/SNP_posns.tsv

# Join matrix of SNPs-by-Affy-ID to the SNP_posns.tsv and report with rows as positions/bases
  join 00_start/SNP_posns.tsv 00_start/SNPs_alleles_w_diploids_arahy3_bySNP.tsv | 
    perl -pe 's/ /\t/g; s/SNP_ID/#SNP_ID/; s/Chr_id/00Chr/' | sort -k2,2 -k3n,3n |
     perl -pe 's/^(\S+)\t(\S+)\t(\S+)\t/$1__$2__$3\t/' > 00_start/SNPs_alleles_rowsSNPs3_by_accn.tsv

  # split by chromosome. Suppress sites missing from duranensis_and_ipaensis, and convert to fasta
  mkdir 00_start/SNPs_alleles_rowsSNPs3_by_accn_CHR
  mkdir 00_start/SNPs_alleles_rowsSNPs3_by_accn_CHR

  for chr in A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 B01 B02 B03 B04 B05 B06 B07 B08 B09 B10; do 
    cat /dev/null > 00_start/SNPs_alleles_rowsSNPs3_by_accn_CHR/SNPs_alleles_rowsSNPs3_$chr.fas
    export CHR=$chr
    cat 00_start/SNPs_alleles_rowsSNPs3_by_accn.tsv |
      perl -lne '$CHR=$ENV{'CHR'}; if ($_=~/#SNP_ID/ || $_=~/$CHR/){print $_}' |
        awk '$1120!~/-/' | transpose.pl | grep -v "#" | perl -pe 's/^(\S+)\t/>$1\n/; s/\t//g' |
        #transpose.pl | grep -v "#" | perl -pe 's/^(\S+)\t/>$1\n/; s/\t//g' |
        cat >> 00_start/SNPs_alleles_rowsSNPs3_by_accn_CHR/SNPs_alleles_rowsSNPs3_$chr.fas
  done

  # Extract ipaensis_duranensis with sequence by chromosome
    grep -A1 ipa 00_start/SNPs_alleles_rowsSNPs3_by_accn_CHR/* | 
      awk '$1!~/^#/ && $1!~/^--/ && $1!~/ipa/' | perl -pe 's/^00.+_(\w\w\w)\.fas-/>$1\n/' > 04_fas/dur_ipa_by_chr.fas

# Join matrix of SNPs-by-Affy-ID to the SNP_posns.tsv and report with rows as accessions
  join 00_start/SNP_posns.tsv 00_start/SNPs_alleles_w_diploids_arahy3_bySNP.tsv | 
    perl -pe 's/ /\t/g; s/SNP_ID/#SNP_ID/; s/Chr_id/00Chr/' | sort -k2,2 -k3n,3n |
     perl -pe 's/^(\S+)\t(\S+)\t(\S+)\t/$1__$2__$3\t/' | transpose.pl > 00_start/SNPs_alleles_rowsSNPs3_by_posn.tsv

# Also make a version where positions are only retained if there are SNPs in the diploids at those positions ... and no scaffolds
  join 00_start/SNP_posns.tsv 00_start/SNPs_alleles_w_diploids_arahy3_bySNP.tsv | 
    perl -pe 's/ /\t/g; s/SNP_ID/#SNP_ID/; s/Chr_id/00Chr/' |
      awk '$2!~/scaff/' | sort -k2,2 -k3n,3n |
      perl -pe 's/^(\S+)\t(\S+)\t(\S+)\t/$1__$2__$3\t/' | 
      awk '$1120!~/-/' |
      transpose.pl > 00_start/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.tsv

# Convert to fasta
  grep -v "^#" 00_start/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.tsv | 
   perl -pe 's/^(\S+)\t/>$1\n/; s/\t//g' > 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.fas

# Compare ipa_dur sequence against all others
  
  makeblastdb -in 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.fas -dbtype nucl -hash_index -out blastdb/SNPs_alleles_rowsSNPs3_by_posn_dip_trim
 
  nohup blastn -query 04_fas/dur_ipa_by_posn.fas \
    -db blastdb/SNPs_alleles_rowsSNPs3_by_posn_dip_trim -evalue 1e-2 -outfmt 6 \
    -num_threads 12 -out blastout/ipa_dur.x.SNPs_alleles_rowsSNPs3_by_posn_dip_trim.bln &


# Do nucmer comparison
  mkdir nucmer_out

  QRY=dur_ipa_by_posn
  REF=SNPs_alleles_rowsSNPs3_by_posn_dip_trim
  
  ~/miniconda2/bin/nucmer --maxmatch -p nucmer_out/$REF.x.$QRY \
    04_fas/$REF.fas 04_fas/$QRY.fas &

  show-coords -rclT nucmer_out/$REF.x.$QRY.delta > nucmer_out/$REF.x.$QRY.coords

  /opt/mummer/3.23/mummerplot -R 04_fas/$REF.fas -Q 04_fas/$QRY.fas --png --vlarge \
    -p nucmer_out/$REF.x.$QRY nucmer_out/$REF.x.$QRY.delta

# Get list of (62) accessions with characteristic identities and coverage vs. duranensis_and_ipaensis
  get_fasta_subset.pl -in 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.fas -lis lis.ipa_dur_matches \
    -out 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_repr.fas

  QRY=dur_ipa_by_chr
  REF=SNPs_alleles_rowsSNPs3_by_posn_dip_trim_repr
  
  ~/miniconda2/bin/nucmer --maxmatch -p nucmer_out/$REF.x.$QRY \
    04_fas/$REF.fas 04_fas/$QRY.fas 

  show-coords -rclT nucmer_out/$REF.x.$QRY.delta > nucmer_out/$REF.x.$QRY.coords

  /opt/mummer/3.23/mummerplot -R 04_fas/$REF.fas -Q 04_fas/$QRY.fas --png --vlarge \
    -p nucmer_out/$REF.x.$QRY nucmer_out/$REF.x.$QRY.delta

# Get fasta of unique accessions, based on list of uniques from peanut_core_v10.xlsx / uniques
  get_fasta_subset.pl -in 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim.fas -lis lis.uniques \
    -out 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.fas

# Make comparisons vs duranensis_ipaensis, and "paint" by similarity
  mkdir 06_haplo_paint
  paint_align_identity.pl -in 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.fas -ref duranensis_and_ipaensis \
  > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.txt

  cat 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.txt | ~/bin/paint_collapse5.pl \
    > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.collapse5.txt

  less -S 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.collapse5.txt | awk 'NR % 10 == 0' | less -S

  cat 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.txt | ~/bin/paint_collapse3.pl \
    > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.collapse3.txt

  less -S 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.collapse3.txt | awk 'NR % 10 == 0' | less -S


# Get fasta at 90% identity, for easier visualization of differences
  perl -pe 's/-/N/g' 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uniques.fas > 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uqNs.fas

  nohup vsearch --cluster_fast 04_fas/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uqNs.fas --id 0.90 --fasta_width 0 \
    --centroids 04_clust/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90.fas \
    --consout 04_clust/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cons90.fas --clusterout_id \
    --uc 04_clust/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cls90.uc &

    # Clusters: 197 Size min 1, max 216, avg 5.0
    # Singletons: 141, 14.2% of seqs, 71.6% of clusters

  # Manually add back duranensis_and_ipaensis, which falls into cluster 21, with 26 members (centroid PI240572_2)

# Highlight haplotypes visually
  paint_align_score.pl -in 04_clust/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90.fas \
    -ref_id duranensis_and_ipaensis > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90_score.txt

  paint_collapse.pl -width 3 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90_score.txt \
    > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90_paint3.txt

  paint_collapse.pl -width 5 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90_score.txt \
    > 06_haplo_paint/SNPs_alleles_rowsSNPs3_by_posn_dip_trim_uq_cen90_paint5.txt

  # See results in haplo90_painted_v02.xlsx

##########
# Investigate SNPs in duranensis and ipaensis again.

# Reshape the BLAST output of the flanking sequences. These had two query sequences per SNP location: one for each variant.
# Write out to file:   diploid_SNP_matches.tsv

  cat blastout/48k_SNPs_combined.x* | perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/; s/gnm2.Aradu/gnm2/' | 
    sort -k1,1 -k4nr,4nr | 
    awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $10, $11 }' 
      AX-147207638  C aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 ACTGTTAAGTGCTTGTGTAGCATATTTTTAAGAATCGGAACATCCAATGACTTTTTGATGGAGCCGAGAAC ACTGTTAAGTGCTTGTGTAGCATATTTTTAAGAATCGGAACATCCAATGACTTTTTGATGGAGCCGAGAAC 231428  231498
      AX-147207665  T aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 CGATCCAACTCTCTCGCGTATTGTTTCATTTCTTCTACCAGAGGAGTCAACGCATCTAAGGTCTCCTTCTT CGATCCAACTCTCTCGCGTATTGTTTCATTTCTTCTACCAGAGGAGTCAACGCATCTAAGGTCTCCTTCTT 413834  413904
      AX-147207665  T araip.K30076.gnm2.B01 100.000 69  0 0 3 71  1.31e-28  128 ATCCAACTCTCTCGCGTATTGTTTCATTTCTTCTACCAGAGGAGTCAACGCATCTAAGGTCTCCTTCTT ATCCAACTCTCTCGCGTATTGTTTCATTTCTTCTACCAGAGGAGTCAACGCATCTAAGGTCTCCTTCTT 17130847  17130779
      AX-147207826  A araip.K30076.gnm2.B01 100.000 66  0 0 6 71  6.09e-27  122 TTTAGCTTGGTCATAAGTGGCTGATGTTTTATCTTTAGCGGTCATAACGGTTTCTCCGATAGCATT  TTTAGCTTGGTCATAAGTGGCTGATGTTTTATCTTTAGCGGTCATAACGGTTTCTCCGATAGCATT  14857501  14857436
      AX-147207826  G aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 TATCTTTTAGCTTGGTCATAAGTGGCTGATGTTTTGTCTTTAGCGGTCATAACGGTTTCTCCGATAGCATT TATCTTTTAGCTTGGTCATAAGTGGCTGATGTTTTGTCTTTAGCGGTCATAACGGTTTCTCCGATAGCATT 1193603 1193673
      AX-147207852  T aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 CATTGTAACCGCCACCGAATTCTTATCCAATGCTGTTGCAACTCCAATCCGCATTGTAAGCTTCTGACCGT CATTGTAACCGCCACCGAATTCTTATCCAATGCTGTTGCAACTCCAATCCGCATTGTAAGCTTCTGACCGT 1354696 1354766
      AX-147207857  A aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 ATTATAGTGGCCCATATGGGTCTAACTACAACGGGATCAATCGGTGCAAAAACTGCTGTTTCACTGGAGGA ATTATAGTGGCCCATATGGGTCTAACTACAACGGGATCAATCGGTGCAAAAACTGCTGTTTCACTGGAGGA 1381799 1381869
      AX-147207858  G aradu.V14167.gnm2.A01 100.000 71  0 0 1 71  7.78e-30  132 TTTTTAAAAACTACACTAGGATTTTGACATACCTGGTGTTCCCCAAAGTGGTCAATGAGGCCTAAAAGGTC TTTTTAAAAACTACACTAGGATTTTGACATACCTGGTGTTCCCCAAAGTGGTCAATGAGGCCTAAAAGGTC 1385261 1385331
      AX-147207858  T araip.K30076.gnm2.B01 100.000 71  0 0 1 71  1.01e-29  132 TTTTTAAAAACTACACTAGGATTTTGACATACCTGTTGTTCCCCAAAGTGGTCAATGAGGCCTAAAAGGTC TTTTTAAAAACTACACTAGGATTTTGACATACCTGTTGTTCCCCAAAGTGGTCAATGAGGCCTAAAAGGTC 14015752  14015682
      AX-147207860  C araip.K30076.gnm2.B01 100.000 71  0 0 1 71  1.01e-29  132 TCATCATCATCTTTGCACATAGGAGGTTCATTTTCCTCCCTTTGCTCATAGCATTCGTTGGACACAGGTTT TCATCATCATCTTTGCACATAGGAGGTTCATTTTCCTCCCTTTGCTCATAGCATTCGTTGGACACAGGTTT 13991450  13991380


# How many matches align only to one subgenome, or to two, or at higher copy? These are filtered first to 100% identity and
# alignment length >=65 (of 71-base queries).
  cat blastout/48k_SNPs_combined.x* | 
    perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/; s/gnm2.Aradu/gnm2/' | 
    sort -k1,1 -k4nr,4nr | 
    awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $10, $11 }' | 
    cut -f1 | uniq -c | sort -n | awk '{print $1}' | uniq -c
      10160 1
        677 2
         41 3
          6 4
          3 5
          1 6
          1 55
          1 101
          1 4275
      # Most of the cases with 2 or 3 matches are subgenome-distinguishing SNPs. 
      # They may be problematic for SNP-calling in a tetraploid, however.
      # The cases with one match should be OK for SNP-calling in a tetraploid, as only one subgenome has a 100% match.
      
      # Examples with two matches at 100%:
      cat blastout/48k_SNPs_combined.x* |
        perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/; s/gnm2.Aradu/gnm2/' | sort -k1,1 -k4nr,4nr | 
        awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {
          print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $10, $11 }' |  
        cut -f1 | uniq -c | sort -n | awk '$1==2' | head
          2 AX-147207665
          2 AX-147207826
          2 AX-147207858
          2 AX-147207860
          2 AX-147207862
          2 AX-147207877
          2 AX-147208865
          2 AX-147209099
          2 AX-147210678
          2 AX-147210969

  # How many of the alignments with two matches have the same base in both locations? These are likely due to subgenome invasions or conversions.
    cat blastout/48k_SNPs_combined.x* | 
      perl -pe 's/(AX-\w+)\.(\w)/$1\t$2/; s/gnm2.Aradu/gnm2/' | sort -k1,1 -k4nr,4nr | 
      awk -v OFS="\t" '$3!~/Scaff/ && $5>=65 && $4==100 {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $10, $11 }' | wc -l
     # 291   -- this is quite a large proportion of all with two matches (677)









