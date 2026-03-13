## *Dopamine signaling enriched striatal gene set predicts striatal dopamine synthesis and physiological activity in vivo

The polygenic architecture of schizophrenia implicates several molecular pathways involved in synaptic function. However, it is unclear how polygenic risk funnels through these pathways to translate into syndromic illness. Using tensor decomposition, we analyze gene co-expression in the caudate nucleus, hippocampus, and dorsolateral prefrontal cortex of post-mortem brain samples from 358 individuals. We identify a set of genes predominantly expressed in the caudate nucleus and associated with both clinical state and genetic risk for schizophrenia that shows dopaminergic selectivity. A higher polygenic risk score for schizophrenia parsed by this set of genes predicts greater dopamine synthesis in the striatum and greater striatal activation during reward anticipation. These results translate dopamine-linked genetic risk variation into in vivo neurochemical and hemodynamic phenotypes in the striatum that have long been implicated in the pathophysiology of schizophrenia.

Sportelli, L., Eisenberg, D.P., Passiatore, R. et al. Dopamine signaling enriched striatal gene set predicts striatal dopamine synthesis and physiological activity in vivo. Nat Commun 15, 3342 (2024). https://doi.org/10.1038/s41467-024-47456-5

---
### Tensor decomposition 
For tensor decomposition analysis of gene co-expression across multiple brain tissues we employed [SDA](https://jmarchini.org/software/#sda).

#### SDA parameters used 
```
sda_static_linux \
	--data /dcs05/lieber/gpergola/Leonardo/caud_dlpfc_hippo_ALL.txt \
	--N 238 \
	--out /dcs05/lieber/gpergola/Leonardo/SDA/SDA.run \
	--max_iter 3000 \
	--num_comps 238 \
	--set_seed 1 2 \
	--num_blocks 20 \
	--num_openmp_threads 20 \
	--remove_zero_comps true
```

To ensure the replicability of the results reported in this work, all processed and aggregated data (Supplementary Data 1–4), SDA input data, and SNPs used to compute C80-PRSs are also available at: https://doi.org/10.5281/zenodo.10699265
