rm collection.sdf
rm -rf data
mkdir data
cd ./data
split -l 100 ../collection.smi
cd ..

for file in $(ls ./data/*)
do
    echo $file
    time xedex -0 -m 1 -i  l -o s $file >> $file'.sdf' &
#    time xedconvert -F $file | xedmin -i x -N -o s | time xedex -o s -m 1 >>  $file'.sdf'  &
done

wait

for file in  $(ls ./data/*.sdf)
do
    echo $file
    cat $file >> collection.sdf
done
