# Multimodal-Plasmacell_manuscript


## Abstract

Antibody-secreting cells are components of humoral immunity by secreting antibodies and providing protection against pathogens. These cells can be of IgM, IgA, or IgG subclass and migrate to class-specific niches. These cells' localisation and rareness make it challenging to define subclass-specific molecular hallmarks. Here, we describe how in-vitro differentiation of peripheral B-cells results in antibody-secreting cells. We use a single-cell multi-modal sequencing approach to find subclass-specific hallmark transcriptional profiles, surface protein expression and signaling pathway activity.

## In this repository

The pages contain code to process, analyze and create figures presented in the full manuscript:

* [QC](QC.html) processes count tables (available at [GSE189953](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189953)): shows how to select high-quality cells; normalize and scale counts.   
* [PC haracterization](cellstate.html) describes the multi-modal characterization of 11-day-old differentiated B-cells to antibody-secreting cells.   
* [Hallmarks Ig-classes](hallmarks.html) describes the identification of Ig-classes and corresponding mRNA, surface protein levels and signaling activity: (phospho-)protein levels.   


## Attribution


We are very thankful for the efforts made by developers of [Seurat](https://satijalab.org/seurat/index.html), [MOFA+](https://biofam.github.io/MOFA2/) and  [workflowr](https://github.com/jdblischak/workflowr). These (well-documented) R-packages enable respectively extensive multi-modal data analysis and reproducible code documentation. 

------ 

The content in this repository is available under the [CC BY 4.0](https://github.com/vanbuggenum/Multimodal-Plasmacell_manuscript/blob/master/LICENSE.md) license. 

For proper attribution, please [cite](https://github.com/vanbuggenum/Multimodal-Plasmacell_manuscript/blob/master/CITATION.bib) our [publication](https://doi.org/10.1016/j.mcpro.2023.100492) containing description and analysis of all presented data and results.


------


