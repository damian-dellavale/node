# NODE

Nested Outlier Detection (NODE) algorithm.

# Description

The NODE algorithm define transient events in time series by detecting anomalies (i.e. outliers) of amplitude across the frequency bands of interest.
For this, the NODE algorithm uses the Local False Discovery Rate (LFDR) method to identify the transient excursions of amplitudes corresponding to outliers,
while controlling for the proportion of false positives.

# Testing the NODE algorithm

1. Download the ```node``` folder in your computer (use the green button ```Code``` above).
2. Download the header ```SubjectUCI29_hdr.mat``` and the data ```SubjectUCI29_data.mat``` files (see the zenodo link below), and copy them into the folder ```node/data/```.
3. Run the ```test_node_v1.m``` script. 
4. Find the results in the folder ```node/matlab/test_node_v1/```.

Note that this is just a toy example showing how to run the NODE algorithm on real iEEG data.
In [Dellavale et al. 2023](https://doi.org/10.1101/2023.04.05.23288085), 
we elaborated on top of these scripts to study the fast-ultradian dynamics of the rate of interictal events observed in iEEG recordings from epileptic patients. 

# References

- Paper describing the NODE algorithm:\
Dellavale, D., Bonini, F., Pizzo, F., et al. Spontaneous fast-ultradian dynamics of polymorphic interictal events in drug-resistant focal epilepsy, Epilepsia (submitted 2023).\
Preprint freely available at DOI: [10.1101/2023.04.05.23288085](https://doi.org/10.1101/2023.04.05.23288085)

- Paper describing the example iEEG dataset:\
Stolk, A., Griffin, S., van der Meij, R. et al. Integrated analysis of anatomical
and electrophysiological human intracranial data. Nat Protoc 13, 1699–1723 (2018).\
DOI: [10.1038/s41596-018-0009-6](https://doi.org/10.1038/s41596-018-0009-6)\
URL: [https://www.nature.com/articles/s41596-018-0009-6](https://www.nature.com/articles/s41596-018-0009-6)

- Link to the iEEG dataset:\
[https://zenodo.org/record/1201560#.ZCxavXbP2Ul](https://zenodo.org/record/1201560#.ZCxavXbP2Ul)

- Link to the Fieldtrip tutorial using the iEEG dataset:\
[https://www.fieldtriptoolbox.org/tutorial/human_ecog/](https://www.fieldtriptoolbox.org/tutorial/human_ecog/)
