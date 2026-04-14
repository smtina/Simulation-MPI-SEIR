import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Configuration des dossiers
DATA_DIR = 'data'
GRAPH_DIR = 'graph'
os.makedirs(GRAPH_DIR, exist_ok=True)

# Utilisation du style 'ggplot' pour se rapprocher du rendu R (fond gris, grille blanche)
plt.style.use('ggplot') 

def plot_seir_curves(filename, title_suffix=""):
    """
    Génère le graphique de l'évolution SEIR à partir d'un fichier CSV.
    """
    filepath = os.path.join(DATA_DIR, filename)
    if not os.path.exists(filepath):
        print(f" Fichier non trouvé : {filepath}. Passage au suivant.")
        return

    try:
        df = pd.read_csv(filepath)
        plt.figure(figsize=(10, 6))
        
        # Tracé des 4 compartiments avec les couleurs standards
        plt.plot(df['step'], df['S'], label='Sains (S)', color='#1f77b4', linewidth=2)
        plt.plot(df['step'], df['E'], label='Exposés (E)', color='#ff7f0e', linewidth=2)
        plt.plot(df['step'], df['I'], label='Infectés (I)', color='#d62728', linewidth=2)
        plt.plot(df['step'], df['R'], label='Récupérés (R)', color='#2ca02c', linewidth=2)
        
        plt.title(f'Évolution Épidémiologique SEIR - {title_suffix} (730 jours)')
        plt.xlabel('Temps (Jours)')
        plt.ylabel('Nombre d\'individus')
        plt.legend(loc='upper right')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        output_name = f"seir_evolution_{title_suffix.lower().replace(' ', '_')}.png"
        plt.savefig(os.path.join(GRAPH_DIR, output_name), dpi=300)
        plt.close()
        print(f" Courbe SEIR générée : {output_name}")
    except Exception as e:
        print(f" Erreur lors du tracé de {filename}: {e}")

def plot_performance_metrics(perf_filename, label="Modèle"):
    """
    Génère les graphiques de Speedup et d'Efficacité.
    """
    filepath = os.path.join(DATA_DIR, perf_filename)
    if not os.path.exists(filepath):
        print(f" Fichier de performance non trouvé : {filepath}")
        return

    try:
        df = pd.read_csv(filepath, names=['Procs', 'Time'])
        df = df.sort_values('Procs')
        
        
        t_seq = df[df['Procs'] == 1]['Time'].iloc[0]
        df['Speedup'] = t_seq / df['Time']
        df['Efficiency'] = df['Speedup'] / df['Procs']
        
        # --- Graphique 1 : Speedup ---
        plt.figure(figsize=(8, 6))
        plt.plot(df['Procs'], df['Speedup'], 'o-', label=f'Speedup {label}', color='blue')
        plt.plot(df['Procs'], df['Procs'], 'k--', label='Speedup idéal (linéaire)')
        plt.title(f'Analyse du Speedup - {label}')
        plt.xlabel('Nombre de Processus (N)')
        plt.ylabel('Speedup S(N)')
        plt.legend()
        plt.savefig(os.path.join(GRAPH_DIR, f'speedup_{label.lower()}.png'), dpi=300)
        
       
        plt.close('all')
        print(f" Metrics HPC ({label}) générées.")
    except Exception as e:
        print(f" Erreur Performance ({perf_filename}): {e}")

if __name__ == "__main__":
    
    # Séquentiel
    plot_seir_curves('seir_sequential.csv', title_suffix="Séquentiel")
    
    # Centralisé
    plot_seir_curves('seir_data_p4_repro.csv', title_suffix="MPI_Central")
      
    # Speedup
    plot_performance_metrics('speedup_allgather.csv', label="Allgather")
