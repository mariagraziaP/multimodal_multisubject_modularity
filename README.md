# multimodal_multisubject_modularity

This repository contains the code to reproduce the multi-modal multi-subject modularity optimization introduced in the paper "Multi-modal and multi-subject modular organization of human brain networks" 

link to the preprint: https://www.biorxiv.org/content/10.1101/2022.01.26.477897v1.abstract



Brief description of the method.

The script contains an extension of the multi-layer modularity optimization, implemented to detect multi-scale modules across subjects and connection modality simultaneously. It optimizes a modularity matrix where on the diagonal blocks we put structural and functional connectiity matrices of a sample of subjects and three resolution parameters {\gamma, \omega, \eta} are present. \gamma is the usual spatial resolution parameter that affects the the dimension and number of communities and can be set so that high \gamma values lead to many small modules and low \gamma values to few big modules. \omega is an inter-layer coupling parameter that links homologous nodes across subjects, while \eta connects nodes across modality. High values of \omega and \eta lead to the detection of communities highly coupled across subjects or connection modality, whereas low values highlight communities that are subject- or modality-specific. The output of the script is a matrix of dimension [Nodes x Subjects x Modality]. 

Multi-scale communities can be found by running the script for different combinations of {\gamma, \omega, \eta}. The set of partitions obtained as an output
