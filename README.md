# Lipaxome-2024

__Auteur__ : Sérine Benachenhou Pilou

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
**Pour vous assurer que Conda est à jour, exécutez :**
```bash
conda update conda
```


## 2- Installation de l'environnement DESeq2

### Création de l'environnement

Si vous voulez utilisez ce pipeline pour utiliser `deseq2` pour l'analyse d'expressions différentielles de gènes, il est compliqué pour Snakemake de créer un environnement virtuel `deseq2` pour faire ce type d'analyse.

Pour y remédier, il faut créer un environnement local où télécharger les packages nécessaires pour lancer deseq2:

```bash
conda create --name deseq2
```

Répondez `y` pour confirmer la location de l'environnement et activez l'environnement deseq2:

```bash
conda activate deseq2
```

### Installation des packages

Ensuite, il faut télécharger les packages suivants dans l'environnement deseq2 (en tenant compte des versions des packages). Si vous avez des problèmes de compatibilités de versions, lisez dans le terminal les packages et les versions proposées pour régler le problèmes.

Il est normal que quand vous lancez une installlation, plusieurs autres packages seront aussi installés. 

**Assurez-vous que les packages installés sont compatibles avec la version de `r-base`.** Vous pouvez trouver toutes les versions d'un package en lançant un `conda search` nom_package.

Voici les packages à installer **dans l'ordre**:

- `r-base` (4.3.2)
```bash
conda install -c conda-forge r-base=4.3.2=h93585b2_0
```

- `r-readr` (2.1.4)
```bash
conda install -c pkgs/r r-readr=2.1.4=r43h884c59f_0
```

- `r-stringr` (1.5.0)
```bash
conda install -c pkgs/r r-stringr=1.5.0=r43h6115d3f_0
```

- `r-xml2` (1.3.5)
```bash
conda install -c conda-forge r-xml2=1.3.5=r43h1ad5fc0_0
```

- `r-tidyverse` (2.0.0)
```bash
conda install -c pkgs/r r-tidyverse=2.0.0=r43h6115d3f_0
```

- `bioconductor-rhdf5` (2.44.0)
```bash
conda install -c bioconda bioconductor-rhdf5=2.44.0=r43hf17093f_1 
```

- `bioconductor-tximport` (1.28.0)
```bash
conda install -c bioconda bioconductor-tximport=1.28.0=r43hdfd78af_0
```

- `yq` (3.2.3)
```bash
conda install -c conda-forge yq=3.2.3=pyhd8ed1ab_0 
```

- `xmltodict` (0.13.0)
```bash
conda install -c pkgs/main xmltodict=0.13.0=py312h06a4308_0
```

- `bioconductor-data-packages` (20230718)
```bash
conda install -c bioconda bioconductor-data-packages=20230718
```

- `bioconductor-tximportdata` (1.28.0)
```bash
conda install -c bioconda bioconductor-tximportdata=1.28.0=r43hdfd78af_0
```

- `r-matrixstats` (1.0.0)
```bash
conda install -c pkgs/r r-matrixstats=1.0.0=r43h76d94ec_0
```

- `bioconductor-matrixgenerics` (1.12.2)
```bash
conda install -c bioconda bioconductor-matrixgenerics=1.12.2=r43hdfd78af_0
```

- `bioconductor-summarizedexperiment` (1.30.2)
```bash
conda install -c bioconda bioconductor-summarizedexperiment=1.30.2=r43hdfd78af_0
```

- `r-lambda.r` (1.2.4)
```bash
conda install -c pkgs/r r-lambda.r=1.2.4=r43h142f84f_0 
```

- `r-futile.logger` (1.4.3)
```bash
conda install -c pkgs/r r-futile.logger=1.4.3=r43h6115d3f_0
```

- `r-bh` (1.81.0_1)
```bash
conda install -c pkgs/r r-bh=1.81.0_1=r43h6115d3f_0  
```

- `r-snow` (0.4_4)
```bash
conda install -c pkgs/r r-snow=0.4_4=r43h142f84f_0  
```

- `bioconductor-biocparallel` (1.34.2)
```bash
conda install -c bioconda bioconductor-biocparallel=1.34.2=r43hf17093f_0   
```

- `r-locfit` (1.5_9.8)
```bash
conda install -c pkgs/r r-locfit=1.5_9.8=r43h76d94ec_0   
```

- `bioconductor-deseq2` (1.40.2)
```bash
conda install -c bioconda bioconductor-deseq2=1.40.2=r43hf17093f_0    
```

### Duplication de l'environnement

Après avoir un environnement deseq2 fonctionnel, vous pouvez dupliquer cet environnement avec un nouveau nom pour y installer snakemake.

```bash
conda create --name smake_deseq2 --clone deseq2
```
Vous devriez avoir au final 3 environnements conda:

- base
- deseq2
- smake_deseq2

## 3- Installation de Mamba

Mamba est un package disponible sur Conda qui permet de télécharger des packages plus rapidement que Conda. Il est conseillé d'utiliser `mamba` pour télécharger `snakemake`.

### Dans l'environnement smake_deseq2

```bash
conda install -c conda-forge mamba
```

### Dans l'environnement base

```bash
conda install -n base -c conda-forge mamba
```

## 4- Installation de Snakemake

### Dans l'environnement smake_deseq2

```bash
# Télécharger snakemake via mamba
mamba install -c conda-forge -c bioconda snakemake
```

### Dans l'environnement base

Il faut créer un nouveau environnement sur Conda et télécharger `snakemake` dans cet environnement. 

```bash
# Télécharger snakemake via mamba
mamba create -c conda-forge -c bioconda -n smake snakemake 
```

### Si vous n'avez pas réussi à télécharger `mamba`:

Vous pouvez aller sur le site de Snakemake pour voir les autres options de téléchargement de `snakemake`:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

## 5- Éxecution du pipeline Snakemake

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
# Lancer un workflow sur un cluster (ou section all du Snakefile)
snakemake --profile ../profile_slurm/

# Lancer la section dge du pipeline sur un cluster
# Assurez-vous d'être dans l'environnement smake_deseq2
snakemake dge --profile ../profile_slurm/

# Lancer la section tgirt du pipeline sur un cluster
snakemake tgirt --profile ../profile_slurm/

# Lancer la section splicing du pipeline sur un cluster
snakemake splicing --profile ../profile_slurm/
```
