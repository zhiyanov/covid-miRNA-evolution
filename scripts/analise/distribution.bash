gisaids=$(cat gisaids.txt)
region="RNA"
for gisaid in $gisaids
do
    echo $gisaid
    cat projects.txt | parallel python3 distribution.py $gisaid $region {}
done
