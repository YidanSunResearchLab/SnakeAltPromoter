


library(rtracklayer)
library(GenomicRanges)

input_gtf  <- "/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation.gtf"  
output_gtf <- "/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/promoters.gtf"
up  <- 50   # upstream bp
down<- 50   # downstream bp


gtf <- import(input_gtf)
tx  <- gtf[gtf$type == "transcript"]


prom <- promoters(tx, upstream = up, downstream = down + 1)

mcols(prom)$gene_id     <- mcols(tx)$gene_id
mcols(prom)$promoter_id <- mcols(tx)$transcript_id

prom$type   <- "promoter"
prom$source <- "HAVANA" 


export(prom, output_gtf, format = "gtf")











featureCounts \
  -a /mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/promoters.gtf \
  -o /mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/cage_promoter_counts.txt \
  -t promoter \
  -g promoter_id \
  -s 1 -M -O \Human_GSM4421327_Human_Left_Ventrice_Tissue_Heart_Failure_CAGE_SINGLE_SRR11351703.Aligned.out.bam






INPUT_BED=/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter.bed
OUTPUT_GTF=/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/cage_promoter_counts.gtf

awk 'BEGIN{OFS="\t"}
{
  chrom=$1; start=$2; end=$3; name=$4; score="."; strand=$6;
  print chrom, "proActiv", "promoter", start+1, end, score, strand, ".", "promoter_id \"" name "\";"
}' "$INPUT_BED" > "$OUTPUT_GTF"
 
 




awk 'BEGIN{OFS="\t"}
     $3=="promoter" {
       start=$4; end=$5;
       if ($7=="+") { new_start=start-100; new_end=start+100 }
       else         { new_start=end-100;   new_end=end+100 }
       if (new_start<1) new_start=1;
       $4=new_start; $5=new_end;
       print
     }' cage_promoter_counts.gtf > cage_promoter_regions.gtf

awk 'BEGIN{OFS="\t"}
     $3=="promoter" {
       start=$4; end=$5;
       if ($7=="+") { new_start=start-100; new_end=start+100 }
       else         { new_start=end-100;   new_end=end+100 }
       if (new_start<1) new_start=1;
       $4=new_start; $5=new_end;
       $9 = "promoter_id \"" substr($9, index($9, "\"")+1) # keep original id
       sub(/[ \t]+$/, "", $9)                               # get rid of space at the end
       if (substr($9, length($9)) != ";") $9 = $9 ";"       # ensure ends with ;
       print
     }' cage_promoter_regions.gtf > cage_promoter_regions_clean.gtf