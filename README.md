# Characterizing dysregulations via cell-cell communications in Alzheimer's brains using single-cell transcriptomes

### Abstract
#### Background
Alzheimer's disease (AD) is a devastating neurodegenerative disorder affecting 44 million people worldwide, leading to cognitive decline, memory loss, and significant impairment in daily functioning. The recent single-cell sequencing technology has revolutionized genetic and genomic resolution by enabling scientists to explore the diversity of gene expression patterns at the finest resolution. Most existing studies have solely focused on molecular perturbations within each cell, but cells live in microenvironments rather than in isolated entities. Here, we leveraged the large-scale and publicly available single-nucleus RNA sequencing in the human prefrontal cortex to investigate cell-to-cell communication in healthy brains and their perturbations in AD. We uniformly processed the snRNA-seq with strict QCs and labeled canonical cell types consistent with the definitions from the BRAIN Initiative Cell Census Network. From ligand and receptor gene expression, we built a high-confidence cell-to-cell communication network to investigate signaling differences between AD and healthy brains.
#### Results
Specifically, we first performed broad communication pattern analyses to highlight that biologically related cell types in normal brains rely on largely overlapping signaling networks and that the AD brain exhibits the irregular inter-mixing of cell types and signaling pathways. Secondly, we performed a more focused cell-type-centric analysis and found that excitatory neurons in AD have significantly increased their communications to inhibitory neurons, while inhibitory neurons and other supporting cells globally decreased theirs to all cells. Then, we delved deeper with a signaling-centric view, showing that canonical signaling pathways CSF, TGFβ, and CX3C are significantly dysregulated in their signaling to the cell type microglia/PVM and from endothelial to neuronal cells for the WNT pathway. Finally, after extracting 23 known AD risk genes, our intracellular communication analysis revealed a strong connection of extracellular ligand genes APP, APOE, and PSEN1 to intracellular AD risk genes TREM2, ABCA1, and APP in the communication from astrocytes and microglia to neurons.
#### Conclusions
In summary, with the novel advances in single-cell sequencing technologies, we show that cellular signaling is regulated in a cell-type-specific manner and that improper regulation of extracellular signaling genes is linked to intracellular risk genes, giving the mechanistic intra- and inter-cellular picture of AD.

### Methodology
1. scRNA-seq Processing:
2. Cell Type Annotation
3. Intercellular Communication Analyses
4. AD Risk Gene Extraction 
5. Intracellular Communication Analyses

<i>For more details, please refer to our research paper.</i>

### Results [Brief]
Results from our C2C analysis have shown that there is a global C2C communication pattern intermixing (inhibitory Chandelier cells in the excitatory group Fig. 2) and signaling pathway misusage in AD brains (e.g. ANGPT, WNT pathways, Fig. 2). Additionally, we also observed a large degree of C2C communication disruption heterogeneity across various cell types. For example, excitatory neurons tend to solely increase their communication strength with inhibitory neurons, while supporting cells and inhibitory neurons globally decrease their communication to most cell types (Fig. 3).  This signifies the importance of employing single-cell technologies in AD studies to dissect the extensive genetic heterogeneity in complex tissues like the human brain. Furthermore, we highlighted the involvement of the neural inflammatory and neural protective pathways, such as WNT, CSF, TGFβ, and CX3C, in AD patients (Fig. 4). Their disturbed behavior can pass erroneous information both inter- and intra-cellularly to directly impact well-known AD risk genes (Fig. 5). 

### System Requirements
Since both conda and R environments were used to install packages for this project, observe that there are two files: `c2c_paper_all_R_pkgs.csv` and `c2c_paper_conda.yml`. To create a new environment, please create it by first using the yml file. Once the env is create, please activate it and then execute the following R code to install all other R packages. 

```
packages <- read.table("path/to/r-packages-list.txt", stringsAsFactors=FALSE)[,1]
install.packages(packages)
```
