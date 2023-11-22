#!/usr/bin/env python

import pandas as pd
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf
import numpy as np
import matplotlib.pyplot as plt

# Load GO terms and associations
obodag = obo_parser.GODag("http://purl.obolibrary.org/obo/go/go-basic.obo")  # Download go-basic.obo from the provided URL
associations = read_gaf("http://current.geneontology.org/annotations/goa_human.gaf.gz")  # Download goa_human.gaf.gz from the provided URL

# Get list of genes
genes = list(pd.read_csv(snakemake.input.genes, sep='\t')['gene'])

# Run GO enrichment analysis
go_enrichment = GOEnrichmentStudy(
    geneids=genes,
    pop=associations.keys(),
    assoc=associations,
    obo=obodag,
    propagate_counts=True,
    alpha=0.05,
    methods=['fdr_bh']
)

# Print significant GO terms
significant_go_terms = [rec.GO for rec in go_enrichment.results if rec.p_fdr_bh < 0.05]

# Plotting GO bar charts
if significant_go_terms:
    # Extract GO terms of interest
    go_terms_of_interest = {term: associations[term] for term in significant_go_terms}

    # Represent GO results with a bar chart
    plt.figure(figsize=(9, 4))

    go_concat = pd.concat(go_terms_of_interest.values(), ignore_index=True)

    bars = pd.DataFrame(
        {'source': go_concat['NS'].values.tolist(),
         'p_val': -np.log10(go_concat['P'])},
        index=go_concat['GO'].values.tolist())

    bars['p_val'].plot(kind="barh", color=bars['source'].replace({'BP': 'indianred', 'MF': 'limegreen', 'CC': 'steelblue'}),
                      width=0.7)
    plt.title('GO Analysis of genes')
    plt.xlabel('$-log_{10}$(p-value)')

    labels = list({'BP': 'indianred', 'MF': 'limegreen', 'CC': 'steelblue'}.keys())
    handles = [plt.Rectangle((0, 0), 1, 1, color=label) for label in labels]
    plt.legend(handles, labels, loc='lower right')

    plt.tight_layout()
    plt.savefig(snakemake.output.bar_chart)
else:
    print("No significant GO terms found.")

