for file in  $(find ./data3/ -name "*.sdf")
do
    echo $file
    cat $file >> extend.sdf
done
