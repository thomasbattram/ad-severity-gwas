# ---------------------------------------------------------------------------
# Simplify GP scripts data
# ---------------------------------------------------------------------------

## This script extracts all GP prescriptions for AD cases available in UK biobank,
## puts them in a table, and removes duplicates. This data will then be passed to clinicians
## so that they can pick out the relevant prescriptions to define severe (or moderate) AD cases

## Probably best to run this script manually on a compute node and tmux to prevent timing out
# interactive-job

## manual args
PHENO=data/ad-ids.txt ## An inclusive set of people with eczema
GPSCRIPTS=data/ukb-pheno/gp_scripts.txt
OUTPREFIX="data/gp_scripts_ecz"

awk -F"\t" 'BEGIN{OFS="\t"}{if(F==1){a[$1]=$1}if(F==2){
 if($1 in a){
  if($4==""){$4="-"}
  if($5==""){$5="-"}
  if($6==""){$6="-"}
  if($7==""){$7="[missing]"}
  print $4,$5,$6,$7}
 }}' F=1 ${PHENO} F=2 ${GPSCRIPTS} \
> ${OUTPREFIX}_inclusive.txt

## output missing data
grep "missing" ${OUTPREFIX}_inclusive.txt > ${OUTPREFIX}_missing_drug.txt

sed 's/ /_/g' ${OUTPREFIX}_inclusive.txt \
> ${OUTPREFIX}_inclusive.nospace.txt

sed 's/[0-9]\+_mg/(!mg_dose!)/g' ${OUTPREFIX}_inclusive.nospace.txt |\
 sed 's/[0-9]\+mg/(!mg_dose!)/g' | sed 's/[0-9]\+_%/(!%_conc!)/g' |\
 sed 's/[0-9]\+%/(!%_conc!)/g' | sed 's/[0-9]\+gms/(!gms_weight!)/g' |\
 sed 's/[0-9]\+ gms/(!gms_weight!)/g' | sed 's/[0-9]\+cm/(!cm_length!)/g' |\
 sed 's/[0-9]\+ cm/(!cm_length!)/g' | sed 's/[0-9]\+ml/(!ml_vol!)/g' |\
 sed 's/[0-9]\+ ml/(!ml_vol!)/g' > ${OUTPREFIX}_inclusive.nodose.txt

sort -k4,4 -k1,3 ${OUTPREFIX}_inclusive.nodose.txt \
 > ${OUTPREFIX}_inclusive.sorted

uniq -c ${OUTPREFIX}_inclusive.sorted |\
 awk 'BEGIN{OFS="\t"}{$1=$1; print $0}' > ${OUTPREFIX}_inclusive.uniq

awk 'BEGIN{n="total_count"; drug="drug"; code2="read_2_concat"; code3="bnf_code_concat"; code4="dmd_code_concat"; OFS="\t"}{
 if($5==drug){
  ## We have a repeat line
  n=n+$1
  code2=code2";"$2
  code3=code3";"$3
  code4=code4";"$4
 } else {
  ## Print out the previous accumulated lines
  print n,code2,code3,code4,drug
  ## Set the new variables to be compared and printed next time
  n=$1
  code2=$2
  code3=$3
  code4=$4
  drug=$5
  drug_concat=substr($5,1,20)
 }
 }END{print n,code2,code3,code4,drug,drug_concat}' ${OUTPREFIX}_inclusive.uniq \
> ${OUTPREFIX}_inclusive.uniq.grouped

awk -F"\t" '{print $0 (NR>1 ? FS substr($5,1,20) : "")}' ${OUTPREFIX}_inclusive.uniq.grouped \
> ${OUTPREFIX}_inclusive.uniq.grouped.concat

rm ${OUTPREFIX}_inclusive.nodose.txt ${OUTPREFIX}_inclusive.nospace.txt \
 ${OUTPREFIX}_inclusive.sorted ${OUTPREFIX}_inclusive.uniq

