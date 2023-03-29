gisaids=$(cat gisaids.txt)
region="RNA"
projects=$(cat projects.txt)

for gisaid in $gisaids
do
    for project in $projects
    do
        # echo $gisaid $region $project
        python3 pvalue.py $gisaid $region $project
    done
done
