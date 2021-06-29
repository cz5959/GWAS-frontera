#!/bin/sh

ID=30770

meta_path=/corral-repl/utexas/Recombining-sex-chro/ukb/data/metadata

field="$ID"; head -n1 $meta_path/ukb45020.txt | tr "\t" "\n" | grep -n -w $field

#test
head -10 $meta_path/ukb45020.txt | awk -F "\t" '{print $13780}'

awk -F "\t" '{print $13781}' $meta_path/ukb45020.txt | paste -d "\t" ids.txt - | awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' - > pheno_RBC_count.txt

