#!/bin/sh

ID=30800

meta_path=/corral-repl/utexas/Recombining-sex-chro/ukb/data/metadata

field="$ID"; head -n1 $meta_path/ukb45020.txt | tr "\t" "\n" | grep -n -w $field

#test
head -10 $meta_path/ukb45020.txt | awk -F "\t" '{print $13780}'

awk -F "\t" '{print $14588}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_protein_total.txt
awk -F "\t" '{print $14338}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_calcium.txt
awk -F "\t" '{print $14324}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_urea.txt
awk -F "\t" '{print $14546}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_SHBG.txt
awk -F "\t" '{print $11666}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_whole_body_fat_mass.txt
awk -F "\t" '{print $8996}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_FVC_best.txt
awk -F "\t" '{print $14436}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_HbA1c.txt


for f in *_h2_Cahoy.results; do mv "$f" "$(echo "f" | sed s/_h2_/_/)"; done
for f in *_h2_Cahoy.log; do mv "$f" "$(echo "f" | sed s/_h2_/_/)"; done