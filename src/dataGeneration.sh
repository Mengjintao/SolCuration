#Compiling cure.cpp
g++ cure.cpp -o cure

mkdir ../newdata
mkdir ../newdata/cure
mkdir ../newdata/clean

#Generate curated PHYS dataset 
./cure ../org/phys/phys_stand.csv
mv collection.csv ../newdata/cure/phys_cure.csv

#Generate curated ESOL dataset 
./cure ../org/phys/phys_stand.csv ../org/aqua/aqua_stand.csv ../org/esol/esol_stand.csv
mv collection.csv ../newdata/cure/esol_cure.csv

#Generate curated AQUA dataset 
./cure ../org/phys/phys_stand.csv ../org/esol/esol_stand.csv ../org/aqua/aqua_stand.csv
mv collection.csv ../newdata/cure/aqua_cure.csv

#Generate curated OCHEM dataset 
./cure ../org/phys/phys_stand.csv ../org/esol/esol_stand.csv ../org/aqua/aqua_stand.csv ../org/ochem/ochem_stand.csv 
mv collection.csv ../newdata/cure/ochem_cure.csv

#Generate curated AQSOL dataset 
./cure ../org/phys/phys_stand.csv ../org/esol/esol_stand.csv ../org/aqua/aqua_stand.csv ../org/ochem/ochem_stand.csv ../org/chembl/chembl_stand.csv ../org/aqsol/aqsol_stand.csv 
mv collection.csv ../newdata/cure/aqsol_cure.csv

#Generate curated CHEMBL dataset 
./cure ../org/phys/phys_stand.csv ../org/esol/esol_stand.csv ../org/aqua/aqua_stand.csv ../org/ochem/ochem_stand.csv ../org/chembl/chembl_stand.csv
mv collection.csv ../newdata/cure/chembl_cure.csv

#Generate curated KINECT dataset 
./cure ../org/kinect/kinect_stand.csv 
mv collection.csv ../newdata/cure/kinect_cure.csv

#Collect all clean datasets 
mv *.csv ../newdata/clean/
