#1
cp ../newdata/clean/aqua_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/aqua_rdnorm.npz

cp ../newdata/clean/phys_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/phys_rdnorm.npz

cp ../newdata/clean/aqsol_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/aqsol_rdnorm.npz

cp ../newdata/clean/chembl_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/chembl_rdnorm.npz

cp ../newdata/clean/esol_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/esol_rdnorm.npz

cp ../newdata/clean/kinect_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/kinect_rdnorm.npz

cp ../newdata/clean/ochem_stand.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/clean/ochem_rdnorm.npz

#2
cp ../newdata/cure/aqua_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/aqua_rdnorm.npz

cp ../newdata/cure/phys_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/phys_rdnorm.npz

cp ../newdata/cure/aqsol_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/aqsol_rdnorm.npz

cp ../newdata/cure/chembl_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/chembl_rdnorm.npz

cp ../newdata/cure/esol_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/esol_rdnorm.npz

cp ../newdata/cure/kinect_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/kinect_rdnorm.npz

cp ../newdata/cure/ochem_cure.csv stand.csv
python ./feature/save_features_test.py
mv hello.npz ../newdata/cure/ochem_rdnorm.npz

rm stand.csv
