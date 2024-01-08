# Lipaxome-2024

__Auteur__ : Sérine Benachenhou

__Adresse Courriel__ :  _<bens3307@usherbrooke.ca>_

Ce pipeline permet d'analyser des données de séquençage à l'ARN (RNA-Seq) à l'aide du gestionnaire Snakemake. Il comporte plusieurs branches d'analyses telles que:

- L'expression différentielle des gènes (DGE)
- Analyse d'enrichissement
- Les évènements de l'épissage alternatif
- L'abondance des ARN du transcriptome
- Biostatistique
- Modèle prédictif

## Versions des packages
- `conda` : 23. 10. 0
- `python` : 3. 12. 0
- `mamba` : 1. 5. 3
- `snakemake` : 7. 32. 4
- `r-base` : 4. 3. 2

## 1- Installation du Miniconda
Miniconda3 a besoin d'être installé (https://docs.conda.io/en/latest/miniconda.html)

Miniconda3 est une distribution minimaliste de Conda qui facilite la gestion des environnements et des paquets logiciels pour entre autre les langages Python et R, permettant ainsi une installation et une maintenance simplifiées des outils et bibliothèques nécessaires.

### Pour les utilisateurs Linux:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
```
Lorsque on vous demande: `Miniconda3 will now be installed into this location:`, vérifier que le chemin proposé vous convient (en général, c'est le cas) et cliquez sur `ENTER` (sinon, spécifiez les chemin que vous voulez avant)

À `Do you wish to update your shell profile to automatically initialize conda?`, répondez `yes`

**Redémarrez le terminal**
```bash
source ~/.bashrc
```

Ça devrait vous donner ceci comme ligne de commande: `(base) serinelo@beluga1:~$` (Environnement `base` du Conda)

La version du Conda proposée sur le site n'est pas nécessairement la version la plus récente.
**Pour vous assurer que Conda est à jour, exécutez :**
```bash
conda update conda
```

## 2- Installation de Mamba

Mamba est un package disponible sur Conda qui permet de télécharger des packages plus rapidement que Conda. Il est conseillé d'utiliser `mamba` pour télécharger `snakemake`.

```bash
conda install -n base -c conda-forge mamba
```

## 3- Installation de Snakemake

Il faut créer un nouveau environnement sur Conda et télécharger `snakemake` dans cet environnement. 

### Si vous avez réussi à télécharger `mamba`:

