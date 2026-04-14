#!/bin/bash

# 1. Préparation des dossiers
mkdir -p data graph bin

echo "--- Benchmarks Épidémiologie (2 ans / 730 jours) ---"

# 2. Test Séquentiel
echo "Running Sequential..."
# Si vous avez codé le nombre de pas en dur (N_STEPS) dans infect_seq.c, 
# assurez-vous de l'avoir modifié à 730 avant de compiler.
./bin/infect_seq 
mv seir_sequential.csv data/ 2>/dev/null

# 3. Test Centralisé (Allgather)
rm -f data/speedup_allgather.csv
for p in 1 2 4 8
do
    echo "Running Centralized with $p processes..."
    # CHANGEMENT : On passe 730 comme premier argument
    /usr/bin/time -f "$p,%e" -a -o data/speedup_allgather.csv mpirun --oversubscribe -np $p ./bin/infect_central 730 1
   
    # Rangement du fichier (le nom généré par le C contient '_repro' si l'arg2 est 1)
    mv seir_data_p${p}_repro.csv data/ 2>/dev/null
done

echo "--- Génération des graphiques ---"
# Vérifiez que plots.py cherche bien les fichiers dans data/
python3 plots.py

echo "Terminé ! Les résultats sur 2 ans sont dans /graph"

