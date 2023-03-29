gisaids=$(cat gisaids.txt)
for gisaid in $gisaids
do
    echo $gisaid
    cat regions.txt | parallel python3 count.py $gisaid {}
done
