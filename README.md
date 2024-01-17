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

**SI VOUS VOULEZ UTILISER CE PIPELINE POUR FAIRE DES ANALYSES D'EXPRESSIONS DIFFÉERENTIELLES, SKIPPER L'ÉTAPE 2-3.**

Mamba est un package disponible sur Conda qui permet de télécharger des packages plus rapidement que Conda. Il est conseillé d'utiliser `mamba` pour télécharger `snakemake`.

```bash
conda install -n base -c conda-forge mamba
```

## 3- Installation de Snakemake

Il faut créer un nouveau environnement sur Conda et télécharger `snakemake` dans cet environnement. 

### Si vous avez réussi à télécharger `mamba`:

```bash
# Télécharger snakemake via mamba
mamba create -c conda-forge -c bioconda -n smake snakemake 
```

### Si vous n'avez pas réussi à télécharger `mamba`:

Vous pouvez aller sur le site de Snakemake pour voir les autres options de téléchargement de `snakemake`:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

## 4- Éxecution du pipeline Snakemake

Le pipeline doit être lancé à partir du répertoire **workflow**. 

Dans les récentes versions de Snakemake, vous pouvez vous créer 2 types de profil pour lancer le pipeline:

- **profile_local**: Si les nœuds du cluster n'ont pas accès à Internet, exécutez en premier lieu les tâches nécessitant Internet localement (i.e all downloads) depuis le répertoire du workflow avec : 

```bash
# Lancer la section all_downloads du pipeline localement
snakemake all_downloads --profile ../profile_local/ 

# Lancer un workflow localement
snakemake --profile ../profile_local/ 
```

- **profile_slurm**: Pour lancer le pipeline sur un cluster Slurm, on peut utiliser la commande suivante pour exécuter (mettre en file d'attente) toutes les tâches à la fois à partir du répertoire de workflow. N'oubliez pas de changer l'adresse d'utilisateur de messagerie pour votre propre adresse e-mail dans cluster.yaml (dans le répertoire profile_slurm).

```bash
# Lancer un workflow sur un cluster
snakemake --profile ../profile_slurm/

# Lancer la section dge du pipeline sur un cluster
snakemake dge --profile ../profile_slurm/

# Lancer la section tgirt du pipeline sur un cluster
snakemake tgirt --profile ../profile_slurm/

# Lancer la section splicing du pipeline sur un cluster
snakemake splicing --profile ../profile_slurm/
```

## 5- Particularités pour la section DGE du pipeline

Pour lancer les tâches **dge** il est plus facile de créer