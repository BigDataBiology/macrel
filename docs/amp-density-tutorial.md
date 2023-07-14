# AMP density tutorial

In this tutorial, we will show how to use the `expected.percontigs` output
produced in the `contigs` and `reads` mode to obtain the density of AMPs
per species.

Consider the expected results in the folder `tests/contigs/expected.percontigs`
as one of the inputs used in this tutorial, and the table for taxonomy of contigs
available in the `example_seqs/example_taxonomy_contigs.tsv.gz`.

First you will need to load tables into your system:

```
import pandas as pd

percontigs = pd.read_table('tests/contigs/expected.percontigs', comment='#')
taxonomy = pd.read_table('example_seqs/example_taxonomy_contigs.tsv.gz')
```

Then, you will need to merge these tables.

```
percontigs = percontigs.merge(on='contig',
                              right=taxonomy,
                              how='outer')
```

Now, we will group results and sum values:

```
percontigs = percontigs.dropna()
percontigs = percontigs.drop('contig', axis=1)
percontigs = percontigs.groupby('taxonomy').agg('sum')
```

By now, you should have a table with the species and the total
of assembled base pairs per species as well as their total number of ORFs,
smORFs and redundant AMPs. Now you just calculate the density as follows:

```
percontigs['AMP_density'] = percontigs.AMPs * 1e6 / percontigs.length
percontigs.to_csv('expected.density',
                  sep='\t',
                  header=True,
                  index=None)
```

The resulting `expected.density` table should be as follows:

| **taxonomy** | **length** | **ORFs** | **smORFs** | **AMPs** | **AMP_density** |
| :---: | :---: | :---: | :---: | :---: | :---: | 
| speciesA | 4208 | 11 | 10 | 0 | 0.000 | 
| speciesB | 11876 | 15 | 8 | 1 | 84.203 |
| speciesC | 5679 | 8 | 6 | 0 | 0.000 |

