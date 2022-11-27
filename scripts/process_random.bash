for i in $(seq 1 600); do
    echo $i
    python3 ../miRNA/rand.py > /dev/null 2>&1
    # python3 process_all.py 2> /dev/null # > /dev/null 2>&1
    # python3 analise_random.py >> ../results/russia_mos_vocs_random_expression.txt 2> /dev/null
    python3 analise_random.py >> ../results/russia_mos_vocs_random.txt 2> /dev/null
done
