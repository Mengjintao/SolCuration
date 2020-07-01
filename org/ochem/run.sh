data=$1
output=standtmp.csv

rm $output
wc -l $data
cat $1 | grep "Train" > input.csv
wc -l input.csv

echo "SMILES,CASRN,EXTERNALID,N,NAME,NAME,ARTICLEID,PUBMEDID,PAGE,TABLE,Water solubility {measured},UNIT {Water solubility},\"Water solubility {measured, converted}\",UNIT {Water solubility},Dataset,Temperature,UNIT {Temperature},Ionic strength,UNIT {Ionic strength},comment (chemical),source,pH,UNIT {pH},Quality code,UNIT {Quality code}" > $output
cat input.csv | grep -v "\[#6\]" | grep -v "\[#7\]" | sort | uniq | grep -v "smiles" >> $output
wc -l $output
python standard_ochem.py
wc -l ochem_stand.csv
#mv ochem_stand.csv ochem_train.csv
rm input.csv $output 
